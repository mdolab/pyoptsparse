# /bin/env python
"""
pyPSQP - the pyPSQP wrapper
"""
# =============================================================================
# PSQP Library
# =============================================================================
try:
    from . import psqp
except ImportError:
    psqp = None
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
# PSQP Optimizer Class
# =============================================================================
class PSQP(Optimizer):
    """
    PSQP Optimizer Class - Inherited from Optimizer Abstract Class
    """

    def __init__(self, raiseError=True, *args, **kwargs):
        name = "PSQP"
        category = "Local Optimizer"
        self.defOpts = {
            "XMAX": [float, 1e16],  # Maximum Stepsize
            "TOLX": [float, 1e-16],  # Variable Change Tolerance
            "TOLC": [float, 1e-6],  # Constraint Violation Tolerance
            "TOLG": [float, 1e-6],  # Lagrangian Gradient Tolerance
            "RPF": [float, 1e-4],  # Penalty Coefficient
            "MIT": [int, 1000],  # Maximum Number of Iterations
            "MFV": [int, 2000],  # Maximum Number of Function Evaluations
            "MET": [int, 2],  # Variable Metric Update (1 - BFGS, 2 - Hoshino)
            "MEC": [int, 2],  # Negative Curvature Correction (1 - None, 2 - Powell's Correction)
            "IPRINT": [int, 2],  # Output Level (0 - None, 1 - Final, 2 - Iter)
            "IOUT": [int, 6],  # Output Unit Number
            "IFILE": [str, "PSQP.out"],  # Output File Name
        }
        self.informs = {
            1: "Change in design variable was less than or equal to tolerance",
            2: "Change in objective function was less than or equal to tolerance",
            3: "Objective function less than or equal to tolerance",
            4: "Maximum constraint value is less than or equal to tolerance",
            11: "Maximum number of iterations exceeded",
            12: "Maximum number of function evaluations exceeded",
            13: "Maximum number of gradient evaluations exceeded",
            -6: "Termination criterion not satisfied, but obtained point is acceptable",
        }
        Optimizer.__init__(self, name, category, self.defOpts, self.informs, *args, **kwargs)

        if psqp is None:
            if raiseError:
                raise Error("There was an error importing the compiled psqp module")

        Optimizer.__init__(self, name, category, self.defOpts, self.informs, *args, **kwargs)

        # PSQP needs Jacobians in dense format
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
            derivatives with finite differenes use 'FD'. 'sens'
            may also be 'CS' which will cause pyOptSpare to compute
            the derivatives using the complex step method. Finally,
            'sens' may be a python function handle which is expected
            to compute the sensitivities directly. For expensive
            function evaluations and/or problems with large numbers of
            design variables this is the preferred method.

        sensStep : float
            Set the step size to use for design variables. Defaults to
            1e-6 when sens is 'FD' and 1e-40j when sens is 'CS'.

        sensMode : str
            Use 'pgc' for parallel gradient computations. Only
            available with mpi4py and each objective evaluation is
            otherwise serial

        storeHistory : str
            File name of the history file into which the history of
            this optimization will be stored

        hotStart : str
            File name of the history file to "replay" for the
            optimziation.  The optimization problem used to generate
            the history file specified in 'hotStart' must be
            **IDENTICAL** to the currently supplied 'optProb'. By
            identical we mean, **EVERY SINGLE PARAMETER MUST BE
            IDENTICAL**. As soon as he requested evaluation point
            from PSQP does not match the history, function and
            gradient evaluations revert back to normal evaluations.

        storeSens : bool
            Flag sepcifying if sensitivities are to be stored in hist.
            This is necessay for hot-starting only.
            """

        self.callCounter = 0
        self.storeSens = storeSens

        if len(optProb.constraints) == 0:
            self.unconstrained = True
            optProb.dummyConstraint = False

        # Set optProb and finalize
        self.optProb = optProb
        self.optProb.finalizeDesignVariables()
        self.optProb.finalizeConstraints()
        self._setInitialCacheValues()
        self._setSens(sens, sensStep, sensMode)
        blx, bux, xs = self._assembleContinuousVariables()
        xi = 3 * np.ones(len(xs), "int")
        xs = np.maximum(xs, blx)
        xs = np.minimum(xs, bux)
        nvar = len(xs)
        ff = self._assembleObjective()

        # pSQP CAN handle two sided constraints, but need to know
        # which ones are which. That is a little tricky to determine,
        # so we will split them into one-sided <= zero constraints and
        # equality=0 constraints such that it is simple to determine
        # the type of constraints. Put the equality constraints frist.

        oneSided = True
        # Set the number of nonlinear constraints snopt *thinks* we have:
        if self.unconstrained:
            ncon = 0
            cf = [0.0]
            cl = []
            cu = []
            ic = []
        else:
            indices, blc, buc, fact = self.optProb.getOrdering(["ne", "le", "ni", "li"], oneSided=oneSided)
            ncon = len(indices)
            cl = np.zeros(ncon)  # -self.getOption('XMAX')*np.ones(ncon)
            cu = np.zeros(ncon)
            cf = np.zeros(ncon + 1)
            ic = 2 * np.ones(ncon, "intc")
            self.optProb.jacIndices = indices
            self.optProb.fact = fact
            self.optProb.offset = buc

            # Also figure out the number of equality and set the type
            # of constraint to 5
            tmp0, __, __, __ = self.optProb.getOrdering(["ne", "le"], oneSided=oneSided)
            ic[0 : len(tmp0)] = 5

        if self.optProb.comm.rank == 0:
            # Set history/hotstart
            self._setHistory(storeHistory, hotStart)

            # ======================================================================
            # PSQP - Objective Values Function
            # ======================================================================
            def pobj(n, x, f):
                if self._checkEval(x):
                    self._internalEval(x)

                f = self.storedData["fobj"]
                psqp.pyflush(self.getOption("IOUT"))
                return f

            # ======================================================================
            # PSQP - Constraint Values Function
            # ======================================================================
            def pcon(n, k, x, g):
                if self._checkEval(x):
                    self._internalEval(x)

                g = self.storedData["fcon"][k - 1]
                psqp.pyflush(self.getOption("IOUT"))
                return g

            # ======================================================================
            # PSQP - Objective Gradients Function
            # ======================================================================
            def pdobj(n, x, df):
                if self._checkEval(x):
                    self._internalEval(x)

                df = self.storedData["gobj"]
                psqp.pyflush(self.getOption("IOUT"))
                return df

            # ======================================================================
            # PSQP - Constraint Gradients Function
            # ======================================================================
            def pdcon(n, k, x, dg):
                if self._checkEval(x):
                    self._internalEval(x)

                dg = self.storedData["gcon"][k - 1]
                psqp.pyflush(self.getOption("IOUT"))
                return dg

            # Setup argument list values
            iterm = np.array(0, int)
            gmax = 0.0
            cmax = 0.0
            opt = self.getOption

            if not opt("IPRINT") <= 2:
                raise Error("Incorrect output level setting (IPRINT) option. Must be <= 2")

            if opt("IPRINT") != 0:
                if os.path.isfile(opt("IFILE")):
                    os.remove(opt("IFILE"))

            # Run PSQP
            t0 = time.time()
            # fmt: off
            psqp.psqp_wrap(nvar, ncon, xs, xi, blx, bux, cf, ic, cl, cu,
                           opt('MIT'), opt('MFV'), opt('MET'), opt('MEC'),
                           opt('XMAX'), opt('TOLX'), opt('TOLC'), opt('TOLG'),
                           opt('RPF'), ff, gmax, cmax, opt('IPRINT'),
                           opt('IOUT'), opt('IFILE'), iterm, pobj, pdobj,
                           pcon, pdcon)
            # fmt: on
            optTime = time.time() - t0

            if opt("IPRINT") > 0:
                psqp.closeunit(opt("IPRINT"))

            # Broadcast a -1 to indcate SLSQP has finished
            self.optProb.comm.bcast(-1, root=0)

            # Store Results
            inform = np.asscalar(iterm)
            sol_inform = {}
            sol_inform["value"] = inform
            sol_inform["text"] = self.informs[inform]
            if self.storeHistory:
                self.metadata["endTime"] = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                self.metadata["optTime"] = optTime
                self.hist.writeData("metadata", self.metadata)
                self.hist.close()

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

    def _on_getOption(self, name):
        pass
