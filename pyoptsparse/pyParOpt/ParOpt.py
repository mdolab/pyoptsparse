import numpy as np
import os
import datetime

try:
    from paropt import ParOpt as _ParOpt
    from mpi4py import MPI
except ImportError:
    _ParOpt = None

from ..pyOpt_optimizer import Optimizer
from ..pyOpt_error import Error


class ParOpt(Optimizer):
    """
    ParOpt optimizer class

    ParOpt has the capability to handle distributed design vectors.
    This is not replicated here since pyOptSparse does not have the
    capability to handle this type of design problem.
    """

    def __init__(self, *args, **kwargs):
        name = "ParOpt"
        category = "Local Optimizer"
        defOpts = {
            "filename": [str, "paropt.out"],
            "algorithm": [str, "ip"],  # Other options: tr
            # Generic options for the interior point method/trust region
            "qn_subspace_size": [int, 10],
            "norm_type": [str, "l2"],  # l1, linfty
            "barrier_strategy": [str, "monotone"],
            "starting_point_strategy": [str, "least_squares_multipliers"],
            "max_iterations": [int, 1000],
            "abs_optimality_tol": [float, 1e-6],
            "rel_function_tol": [float, 0.0],
            "penalty_gamma": [float, 1000.0],
            "barrier_fraction": [float, 0.25],
            "barrier_power": [float, 1.0],
            "reset_hessian_frequency": [int, 100000],
            "bfgs_update_type": [str, "skip"],  # or 'damped'
            "affine_step_multiplier_min": [float, 1.0],
            "init_barrier_parameter": [float, 0.1],
            "max_linesearch_iters": [int, 10],
            "armijo_parameter": [float, 1e-3],
            "penalty_descent_fraction": [float, 0.3],
            "min_penalty_parameter": [float, 0.0],
            "output_level": [int, 0],
            # Trust region specifics
            "tr_init_size": [float, 0.01],
            "tr_max_size": [float, 1.0],
            "tr_min_size": [float, 0.0],
            "tr_eta": [float, 0.25],
            "tr_penalty_gamma": [float, 10.0],
            "tr_max_iterations": [int, 250],
            "tr_abs_optimality_tol": [float, 1e-6],
        }
        informs = {}
        if _ParOpt is None:
            raise Error("There was an error importing ParOpt")

        self.set_options = []
        Optimizer.__init__(self, name, category, defOpts, informs, *args, **kwargs)

        # ParOpt requires a dense Jacobian format
        self.jacType = "dense2d"

        return

    def __call__(
        self, optProb, sens=None, sensStep=None, sensMode=None, storeHistory=None, hotStart=None, storeSens=True
    ):
        """
        This is the main routine used to solve the optimization
        problem.

        Parameters
        ----------
        optProb : Optimization or Solution class instance
            This is the complete description of the optimization problem
            to be solved by the optimizer

        sens : str or python Function.
            Specifiy method to compute sensitivities. To
            explictly use pyOptSparse gradient class to do the
            derivatives with finite differenes use \'FD\'. \'sens\'
            may also be \'CS\' which will cause pyOptSpare to compute
            the derivatives using the complex step method. Finally,
            \'sens\' may be a python function handle which is expected
            to compute the sensitivities directly. For expensive
            function evaluations and/or problems with large numbers of
            design variables this is the preferred method.

        sensStep : float
            Set the step size to use for design variables. Defaults to
            1e-6 when sens is \'FD\' and 1e-40j when sens is \'CS\'.

        sensMode : str
            Use \'pgc\' for parallel gradient computations. Only
            available with mpi4py and each objective evaluation is
            otherwise serial

        storeHistory : str
            File name of the history file into which the history of
            this optimization will be stored

        hotStart : str
            File name of the history file to "replay" for the
            optimziation.  The optimization problem used to generate
            the history file specified in \'hotStart\' must be
            **IDENTICAL** to the currently supplied \'optProb\'. By
            identical we mean, **EVERY SINGLE PARAMETER MUST BE
            IDENTICAL**. As soon as he requested evaluation point
            from ParOpt does not match the history, function and
            gradient evaluations revert back to normal evaluations.

        storeSens : bool
            Flag sepcifying if sensitivities are to be stored in hist.
            This is necessay for hot-starting only.
            """

        self.callCounter = 0
        self.storeSens = storeSens

        if len(optProb.constraints) == 0:
            # If the problem is unconstrained, add a dummy constraint.
            self.unconstrained = True
            optProb.dummyConstraint = True

        # Save the optimization problem and finalize constraint
        # jacobian, in general can only do on root proc
        self.optProb = optProb
        self.optProb.finalizeDesignVariables()
        self.optProb.finalizeConstraints()
        self._setInitialCacheValues()
        self._setSens(sens, sensStep, sensMode)
        blx, bux, xs = self._assembleContinuousVariables()
        xs = np.maximum(xs, blx)
        xs = np.minimum(xs, bux)

        # The number of design variables
        n = len(xs)

        oneSided = True

        if self.unconstrained:
            m = 0
        else:
            indices, blc, buc, fact = self.optProb.getOrdering(["ne", "le", "ni", "li"], oneSided=oneSided)
            m = len(indices)
            self.optProb.jacIndices = indices
            self.optProb.fact = fact
            self.optProb.offset = buc

        if self.optProb.comm.rank == 0:
            # Set history/hotstart
            self._setHistory(storeHistory, hotStart)

            class Problem(_ParOpt.Problem):
                def __init__(self, ptr, n, m, xs, blx, bux):
                    super(Problem, self).__init__(MPI.COMM_SELF, n, m)
                    self.ptr = ptr
                    self.n = n
                    self.m = m
                    self.xs = xs
                    self.blx = blx
                    self.bux = bux
                    self.fobj = 0.0
                    return

                def getVarsAndBounds(self, x, lb, ub):
                    """Get the variable values and bounds"""
                    # Find the average distance between lower and upper bound
                    bound_sum = 0.0
                    for i in range(len(x)):
                        if self.blx[i] <= -1e20 or self.bux[i] >= 1e20:
                            bound_sum += 1.0
                        else:
                            bound_sum += self.bux[i] - self.blx[i]
                    bound_sum = bound_sum / len(x)

                    for i in range(len(x)):
                        x[i] = self.xs[i]
                        lb[i] = self.blx[i]
                        ub[i] = self.bux[i]
                        if self.xs[i] <= self.blx[i]:
                            x[i] = self.blx[i] + 0.5 * np.min((bound_sum, self.bux[i] - self.blx[i]))
                        elif self.xs[i] >= self.bux[i]:
                            x[i] = self.bux[i] - 0.5 * np.min((bound_sum, self.bux[i] - self.blx[i]))

                    return

                def evalObjCon(self, x):
                    """Evaluate the objective and constraint values"""
                    fobj, fcon, fail = self.ptr._masterFunc(x[:], ["fobj", "fcon"])
                    self.fobj = fobj
                    return fail, fobj, -fcon

                def evalObjConGradient(self, x, g, A):
                    gobj, gcon, fail = self.ptr._masterFunc(x[:], ["gobj", "gcon"])
                    g[:] = gobj[:]
                    for i in range(self.m):
                        A[i][:] = -gcon[i][:]
                    return fail

            # Create the ParOpt problem class
            problem = Problem(self, n, m, xs, blx, bux)

            # Get the algorithm/subspace size parameters
            algorithm = self.getOption("algorithm").lower()
            qn_subspace_size = self.getOption("qn_subspace_size")
            filename = self.getOption("filename")

            optTime = MPI.Wtime()
            if algorithm == "ip":
                # Create the optimizer
                opt = _ParOpt.InteriorPoint(problem, qn_subspace_size, _ParOpt.BFGS)

                # Set the ParOpt options
                self._set_paropt_options(opt)

                # Optimize!
                opt.setOutputFile(filename)
                opt.optimize()
            else:
                # Optimality tolerance
                opt_tol = self.getOption("abs_optimality_tol")

                # Trust region algorithm options
                tr_init_size = self.getOption("tr_init_size")
                tr_max_size = self.getOption("tr_max_size")
                tr_min_size = self.getOption("tr_min_size")
                tr_eta = self.getOption("tr_eta")
                tr_penalty_gamma = self.getOption("tr_penalty_gamma")
                tr_opt_abs_tol = self.getOption("tr_abs_optimality_tol")
                tr_max_iterations = self.getOption("tr_max_iterations")

                # Create the quasi-Newton Hessian approximation
                qn = _ParOpt.LBFGS(problem, subspace=qn_subspace_size)
                subproblem = _ParOpt.QuadraticSubproblem(problem, qn)

                # Create the trust region problem
                tr = _ParOpt.TrustRegion(subproblem, tr_init_size, tr_min_size, tr_max_size, tr_eta, tr_penalty_gamma)

                # Create the ParOpt problem
                opt = _ParOpt.InteriorPoint(subproblem, qn_subspace_size, _ParOpt.NO_HESSIAN_APPROX)

                # Set the ParOpt options
                self._set_paropt_options(opt)

                # Set the output file name
                opt.setOutputFile(filename)
                tr.setOutputFile(os.path.splitext(filename)[0] + ".tr")

                # Use the adaptive penalty update scheme by default
                tr.setAdaptiveGammaUpdate(1)
                tr.setPenaltyGammaMax(1e3)

                # Set parameters for the trust-region algorithm
                tr.setMaxTrustRegionIterations(tr_max_iterations)

                # Set the tolerance
                tr.setTrustRegionTolerances(opt_tol, opt_tol, opt_tol)

                # Set optimality tolerance for the trust region problem
                opt.setAbsOptimalityTol(tr_opt_abs_tol)

                # Optimize the problem
                tr.optimize(opt)

                # Get the optimized point
                x, z, zw, zl, zu = opt.getOptimizedPoint()

            # Set the total opt time
            optTime = MPI.Wtime() - optTime

            # Get the obective function value
            fobj = problem.fobj

            # Get the optimized point
            x, z, zw, zl, zu = opt.getOptimizedPoint()

            if self.storeHistory:
                self.metadata["endTime"] = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                self.metadata["optTime"] = optTime
                self.hist.writeData("metadata", self.metadata)
                self.hist.close()

            # Create the optimization solution. Note that the signs on the multipliers
            # are switch since ParOpt uses a formulation with c(x) >= 0, while pyOpt
            # uses g(x) = -c(x) <= 0. Therefore the multipliers are reversed.
            sol_inform = {}

            # If number of constraints is zero, ParOpt returns z as None.
            # Thus if there is no constraints, should pass an empty list
            # to multipliers instead of z.
            if z is not None:
                sol = self._createSolution(optTime, sol_inform, fobj, x[:], multipliers=-z)
            else:
                sol = self._createSolution(optTime, sol_inform, fobj, x[:], multipliers=[])

            # Indicate solution finished
            self.optProb.comm.bcast(-1, root=0)
        else:  # We are not on the root process so go into waiting loop:
            self._waitLoop()
            sol = None

        # Communication solution and return
        sol = self._communicateSolution(sol)

        return sol

    def _set_paropt_options(self, opt):
        """
        set all of the the options in self.set_options in the ipopt instance nlp
        """
        # Set Options from the local options dictionary
        # ---------------------------------------------

        for key in self.options:
            if key != "defaults":
                value = self.getOption(key)

                if key == "norm_type":
                    if value == "l1":
                        opt.setNormType(_ParOpt.L1_NORM)
                    elif value == "linfty":
                        opt.setNormType(_ParOpt.INFTY_NORM)
                    elif value == "l2":
                        opt.setNormType(_ParOpt.L2_NORM)
                elif key == "barrier_strategy":
                    if value == "monotone":
                        opt.setBarrierStrategy(_ParOpt.MONOTONE)
                    elif value == "mehrotra":
                        opt.setBarrierStrategy(_ParOpt.MEHROTRA)
                    elif value == "complementarity_fraction":
                        opt.setBarrierStrategy(_ParOpt.COMPLEMENTARITY_FRACTION)
                elif key == "starting_point_strategy":
                    if value == "none":
                        opt.setStartingPointStrategy(_ParOpt.NO_START_STRATEGY)
                    elif value == "least_squares_multipliers":
                        opt.setStartingPointStrategy(_ParOpt.LEAST_SQUARES_MULTIPLIERS)
                    elif value == "affine_step":
                        opt.setStartingPointStrategy(_ParOpt.AFFINE_STEP)
                elif key == "max_iterations":
                    opt.setMaxMajorIterations(value)
                elif key == "abs_optimality_tol":
                    opt.setAbsOptimalityTol(value)
                elif key == "rel_function_tol":
                    opt.setRelFunctionTol(value)
                elif key == "penalty_gamma":
                    opt.setPenaltyGamma(value)
                elif key == "barrier_power":
                    opt.setBarrierFraction(value)
                elif key == "barrier_power":
                    opt.setBarrierPower(value)
                elif key == "reset_hessian_frequency":
                    opt.setHessianResetFreq(value)
                elif key == "bfgs_update_type":
                    if value == "skip":
                        opt.setBFGSUpdateType(_ParOpt.SKIP_NEGATIVE_CURVATURE)
                    elif value == "damped":
                        opt.setBFGSUpdateType(_ParOpt.DAMPED_UPDATE)
                elif key == "affine_step_multiplier_min":
                    opt.setStartAffineStepMultiplierMin(value)
                elif key == "init_barrier_parameter":
                    opt.setInitBarrierParameter(value)
                elif key == "max_linesearch_iters":
                    opt.setMaxLineSearchIters(value)
                elif key == "armijo_parameter":
                    opt.setArmijoParam(value)
                elif key == "penalty_descent_fraction":
                    opt.setPenaltyDescentFraction(value)
                elif key == "min_penalty_parameter":
                    opt.setPenaltyDescentFraction(value)
                elif key == "output_level":
                    opt.setOutputLevel(value)

    def _on_setOption(self, name, value):
        pass

    def _on_getOption(self, name, value):
        pass
