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
        if _ParOpt is None:
            raise Error("There was an error importing ParOpt")

        # Create and fill-in the dictionary of default option values
        defOpts = {}
        options = _ParOpt.getOptionsInfo()
        for name in options:
            # Get the type and default value of the named argument
            _type = None
            if options[name].option_type == 'bool':
                _type = bool
            elif options[name].option_type == 'int':
                _type = int
            elif options[name].option_type == 'float':
                _type = float
            else:
                _type = str
            default_value = options[name].default

            # Set the entry into the dictionary
            defOpts[name] = [_type, default_value]

        self.set_options = {}
        informs = {}
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
                    """Evaluate the objective and constraint gradients"""
                    gobj, gcon, fail = self.ptr._masterFunc(x[:], ["gobj", "gcon"])
                    g[:] = gobj[:]
                    for i in range(self.m):
                        A[i][:] = -gcon[i][:]
                    return fail

            optTime = MPI.Wtime()

            # Optimize the problem
            problem = Problem(self, n, m, xs, blx, bux)
            optimizer = _ParOpt.Optimizer(problem, self.set_options)
            optimizer.optimize()
            x, z, zw, zl, zu = optimizer.getOptimizedPoint()

            # Set the total opt time
            optTime = MPI.Wtime() - optTime

            # Get the obective function value
            fobj = problem.fobj

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
                sol = self._createSolution(optTime, sol_inform, fobj, x[:],
                                           multipliers=-z)
            else:
                sol = self._createSolution(optTime, sol_inform, fobj, x[:],
                                           multipliers=[])

            # Indicate solution finished
            self.optProb.comm.bcast(-1, root=0)
        else:  # We are not on the root process so go into waiting loop:
            self._waitLoop()
            sol = None

        # Communication solution and return
        sol = self._communicateSolution(sol)

        return sol


    def _on_setOption(self, name, value):
        """
        Add the value to the set_options dictionary.
        """
        self.set_options[name] = value
