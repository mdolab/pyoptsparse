"""
pySLSQP - A variation of the pySLSQP wrapper specificially designed to
work with sparse optimization problems.
"""

# Standard Python modules
import datetime
import os
import time

# External modules
import numpy as np

# Local modules
from ..pyOpt_error import pyOptSparseWarning
from ..pyOpt_optimizer import Optimizer
from ..pyOpt_utils import try_import_compiled_module_from_path

# import the compiled module
THIS_DIR = os.path.dirname(os.path.abspath(__file__))
slsqp = try_import_compiled_module_from_path("slsqp", THIS_DIR, raise_warning=True)


class SLSQP(Optimizer):
    """
    SLSQP Optimizer Class - Inherited from Optimizer Abstract Class
    """

    def __init__(self, raiseError=True, options={}):
        name = "SLSQP"
        category = "Local Optimizer"
        defOpts = self._getDefaultOptions()
        informs = self._getInforms()
        if isinstance(slsqp, str) and raiseError:
            raise ImportError(slsqp)

        self.set_options = []
        super().__init__(name, category, defaultOptions=defOpts, informs=informs, options=options)

        # SLSQP needs Jacobians in dense format
        self.jacType = "dense2d"

    @staticmethod
    def _getDefaultOptions():
        defOpts = {
            "ACC": [float, 1e-6],
            "MAXIT": [int, 500],
            "IPRINT": [int, 1],
            "IOUT": [int, 60],
            "IFILE": [str, "SLSQP.out"],
        }
        return defOpts

    @staticmethod
    def _getInforms():
        informs = {
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
        return informs

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
        self.startTime = time.time()
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
        self.optProb.finalize()
        # Set history/hotstart
        self._setHistory(storeHistory, hotStart)
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
            # =================================================================
            # SLSQP - Objective/Constraint Values Function
            # =================================================================
            def slfunc(m, me, la, n, f, g, x):
                if (x < blx).any() or (x > bux).any():
                    pyOptSparseWarning("Values in x were outside bounds during" " a minimize step, clipping to bounds")
                fobj, fcon, fail = self._masterFunc(np.clip(x, blx, bux), ["fobj", "fcon"])
                f = fobj
                g[0:m] = -fcon
                slsqp.pyflush(self.getOption("IOUT"))
                return f, g

            # =================================================================
            # SLSQP - Objective/Constraint Gradients Function
            # =================================================================
            def slgrad(m, me, la, n, f, g, df, dg, x):
                gobj, gcon, fail = self._masterFunc(np.clip(x, blx, bux), ["gobj", "gcon"])
                df[0:n] = gobj.copy()
                dg[0:m, 0:n] = -gcon.copy()
                slsqp.pyflush(self.getOption("IOUT"))
                return df, dg

            # Setup argument list values
            la = max(m, 1)
            gg = np.zeros([la], float)
            df = np.zeros([n + 1], float)
            dg = np.zeros([la, n + 1], float)
            acc = np.array(self.getOption("ACC"), float)
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
            lw = np.array(lwM, int)
            w = np.zeros(lw, float)
            ljwM = max(mineq, (n + 1) - meq)
            ljw = np.array(ljwM, int)
            jw = np.zeros(ljw, np.intc)
            nfunc = np.array(0, int)
            ngrad = np.array(0, int)

            # Run SLSQP
            t0 = time.time()
            # fmt: off
            slsqp.slsqp(m, meq, la, n, xs, blx, bux, ff, gg, df, dg, acc, maxit,
                        iprint, iout, ifile, mode, w, lw, jw, ljw, nfunc,
                        ngrad, slfunc, slgrad)
            # fmt: on
            optTime = time.time() - t0

            # Clip final result to user bounds (this occurs during the optimization as well
            # so this just makes the output consistent with what the optimizer sees)
            xs = np.clip(xs, blx, bux)

            # some entries of W include the lagrange multipliers
            # for each constraint, there are two entries (lower, upper).
            # if only one is active, look for the nonzero. If both are active, take the first one
            # FIXME: this does not currently work, so we do not save lambdaStar
            # to the solution object
            lambdaStar = []
            idx = 0

            for c_name in optProb.constraints:
                c = optProb.constraints[c_name]
                for _j in range(c.ncon):
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
            inform = mode.item()
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
