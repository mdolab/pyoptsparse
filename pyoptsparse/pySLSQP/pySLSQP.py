# /bin/env python
"""
pySLSQP - A variation of the pySLSQP wrapper specificially designed to
work with sparse optimization problems.
"""
# =============================================================================
# SLSQP Library
# =============================================================================
try:
    from . import slsqp
except ImportError:
    slsqp = None
# =============================================================================
# Standard Python modules
# =============================================================================
import os
import time
import datetime

# =============================================================================
# External Python modules
# =============================================================================
import numpy as np

# ===========================================================================
# Extension modules
# ===========================================================================
from ..pyOpt_optimizer import Optimizer
from ..pyOpt_error import Error

# =============================================================================
# SLSQP Optimizer Class
# =============================================================================
class SLSQP(Optimizer):
    """
    SLSQP Optimizer Class - Inherited from Optimizer Abstract Class
    """

    def __init__(self, raiseError=True, *args, **kwargs):
        name = "SLSQP"
        category = "Local Optimizer"
        self.defOpts = {
            # SLSQP Options
            "ACC": [float, 1e-6],  # Convergence Accurancy
            "MAXIT": [int, 500],  # Maximum Iterations
            "IPRINT": [int, 1],  # Output Level (<0 - None, 0 - Screen, 1 - File)
            "IOUT": [int, 6],  # Output Unit Number
            "IFILE": [str, "SLSQP.out"],  # Output File Name
        }
        self.informs = {
            -1: "Gradient evaluation required (g & a)",
            0: "Optimization terminated successfully.",
            1: "Function evaluation required (f & c)",
            2: "More equality constraints than independent variables",
            3: "More than 3*n iterations in LSQ subproblem",
            4: "Inequality constraints incompatible",
            5: "Singular matrix E in LSQ subproblem",
            6: "Singular matrix C in LSQ subproblem",
            7: "Rank-deficient equality constraint subproblem HFTI",
            8: "Positive directional derivative for linesearch",
            9: "Iteration limit exceeded",
        }
        if slsqp is None:
            if raiseError:
                raise Error("There was an error importing the compiled slsqp module")

        self.set_options = []
        Optimizer.__init__(self, name, category, self.defOpts, self.informs, *args, **kwargs)

        # SLSQP needs Jacobians in dense format
        self.jacType = "dense2d"

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
            from SLSQP does not match the history, function and
            gradient evaluations revert back to normal evaluations.

        storeSens : bool
            Flag sepcifying if sensitivities are to be stored in hist.
            This is necessay for hot-starting only.
            """

        self.callCounter = 0
        self.storeSens = storeSens

        if len(optProb.constraints) == 0:
            # If the user *actually* has an unconstrained problem,
            # slsqp sort of chokes with that....it has to have at
            # least one constraint. So we will add one
            # automatically here:
            self.unconstrained = True
            optProb.dummyConstraint = True

        # Save the optimization problem and finalize constraint
        # Jacobian, in general can only do on root proc
        self.optProb = optProb
        self.optProb.finalizeDesignVariables()
        self.optProb.finalizeConstraints()
        self._setInitialCacheValues()
        self._setSens(sens, sensStep, sensMode)
        blx, bux, xs = self._assembleContinuousVariables()
        xs = np.maximum(xs, blx)
        xs = np.minimum(xs, bux)
        n = len(xs)
        ff = self._assembleObjective()

        oneSided = True
        # Set the number of nonlinear constraints snopt *thinks* we have:
        if self.unconstrained:
            m = 0
            meq = 0
        else:
            indices, blc, buc, fact = self.optProb.getOrdering(["ne", "le", "ni", "li"], oneSided=oneSided)
            m = len(indices)

            self.optProb.jacIndices = indices
            self.optProb.fact = fact
            self.optProb.offset = buc

            # Also figure out the number of equality:
            tmp0, __, __, __ = self.optProb.getOrdering(["ne", "le"], oneSided=oneSided)
            meq = len(tmp0)

        if self.optProb.comm.rank == 0:
            # Set history/hotstart
            self._setHistory(storeHistory, hotStart)

            # =================================================================
            # SLSQP - Objective/Constraint Values Function
            # =================================================================
            def slfunc(m, me, la, n, f, g, x):
                fobj, fcon, fail = self._masterFunc(x, ["fobj", "fcon"])
                f = fobj
                g[0:m] = -fcon
                slsqp.pyflush(self.getOption("IOUT"))
                return f, g

            # =================================================================
            # SLSQP - Objective/Constraint Gradients Function
            # =================================================================
            def slgrad(m, me, la, n, f, g, df, dg, x):
                gobj, gcon, fail = self._masterFunc(x, ["gobj", "gcon"])
                df[0:n] = gobj.copy()
                dg[0:m, 0:n] = -gcon.copy()
                slsqp.pyflush(self.getOption("IOUT"))
                return df, dg

            # Setup argument list values
            la = max(m, 1)
            gg = np.zeros([la], np.float)
            df = np.zeros([n + 1], np.float)
            dg = np.zeros([la, n + 1], np.float)
            acc = np.array([self.getOption("ACC")], np.float)
            maxit = self.getOption("MAXIT")
            iprint = self.getOption("IPRINT")
            iout = self.getOption("IOUT")
            ifile = self.getOption("IFILE")
            if iprint >= 0:
                if os.path.isfile(ifile):
                    os.remove(ifile)

            mode = np.array(0, int)
            mineq = m - meq + 2 * (n + 1)
            lsq = (n + 1) * ((n + 1) + 1) + meq * ((n + 1) + 1) + mineq * ((n + 1) + 1)
            lsi = ((n + 1) - meq + 1) * (mineq + 2) + 2 * mineq
            lsei = ((n + 1) + mineq) * ((n + 1) - meq) + 2 * meq + (n + 1)
            slsqpb = (n + 1) * (n / 2) + 2 * m + 3 * n + 3 * (n + 1) + 1
            lwM = lsq + lsi + lsei + slsqpb + n + m
            lw = np.array([lwM], np.int)
            w = np.zeros(lw, np.float)
            ljwM = max(mineq, (n + 1) - meq)
            ljw = np.array([ljwM], np.int)
            jw = np.zeros(ljw, np.intc)
            nfunc = np.array([0], np.int)
            ngrad = np.array([0], np.int)

            # Run SLSQP
            t0 = time.time()
            # fmt: off
            slsqp.slsqp(m, meq, la, n, xs, blx, bux, ff, gg, df, dg, acc, maxit,
                        iprint, iout, ifile, mode, w, lw, jw, ljw, nfunc,
                        ngrad, slfunc, slgrad)
            # fmt: on
            optTime = time.time() - t0

            # some entries of W include the lagrange multipliers
            # for each constraint, there are two entries (lower, upper).
            # if only one is active, look for the nonzero. If both are active, take the first one
            # FIXME: this does not currently work, so we do not save lambdaStar
            # to the solution object
            lambdaStar = []
            idx = 0

            for c_name in optProb.constraints:
                c = optProb.constraints[c_name]
                for j in range(c.ncon):
                    lambdaStar_lower = w[2 * idx]
                    lambdaStar_upper = w[2 * idx + 1]
                    if abs(lambdaStar_lower) > 1e-100:
                        lambdaStar.append(lambdaStar_lower)
                    else:
                        lambdaStar.append(lambdaStar_upper)
                    idx += 1

            if self.storeHistory:
                self.metadata["endTime"] = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                self.metadata["optTime"] = optTime
                self.hist.writeData("metadata", self.metadata)
                self.hist.close()

            if iprint > 0:
                slsqp.closeunit(self.getOption("IOUT"))

            # Broadcast a -1 to indcate SLSQP has finished
            self.optProb.comm.bcast(-1, root=0)

            # Store Results
            inform = np.asscalar(mode)
            sol_inform = {}
            sol_inform["value"] = inform
            sol_inform["text"] = self.informs[inform]

            # Create the optimization solution
            sol = self._createSolution(optTime, sol_inform, ff, xs)

        else:  # We are not on the root process so go into waiting loop:
            self._waitLoop()
            sol = None

        # Communication solution and return
        sol = self._communicateSolution(sol)

        return sol

    def _on_setOption(self, name, value):
        pass

    def _on_getOption(self, name, value):
        pass
