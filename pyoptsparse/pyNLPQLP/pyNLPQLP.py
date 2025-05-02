"""
pyNLPQLP - A pyOptSparse wrapper for Schittkowski's NLPQLP
optimization algorithm.
"""

# Standard Python modules
import datetime
import os
import time

# External modules
import numpy as np

# Local modules
from ..pyOpt_optimizer import Optimizer
from ..pyOpt_utils import try_import_compiled_module_from_path

# import the compiled module
THIS_DIR = os.path.dirname(os.path.abspath(__file__))
nlpqlp = try_import_compiled_module_from_path("nlpqlp", THIS_DIR)


class NLPQLP(Optimizer):
    """
    NLPQL Optimizer Class - Inherited from Optimizer Abstract Class
    """

    def __init__(self, raiseError=True, options={}):
        name = "NLPQLP"
        category = "Local Optimizer"
        defOpts = self._getDefaultOptions()
        informs = self._getInforms()
        if isinstance(nlpqlp, str) and raiseError:
            raise ImportError(nlpqlp)

        super().__init__(name, category, defaultOptions=defOpts, informs=informs, options=options)
        # NLPQLP needs Jacobians in dense format
        self.jacType = "dense2d"

    @staticmethod
    def _getInforms():
        informs = {
            -2: (
                "Compute gradient values w.r.t. the variables stored in"
                + " first column of X, and store them in DF and DG."
                + " Only derivatives for active constraints ACTIVE(J)=.TRUE. need to be computed."
            ),
            -1: (
                "Compute objective fn and all constraint values subject"
                + "the variables found in the first L columns of X, and store them in F and G."
            ),
            0: "The optimality conditions are satisfied.",
            1: "The algorithm has been stopped after MAXIT iterations.",
            2: "The algorithm computed an uphill search direction.",
            3: "Underflow occurred when determining a new approximation matrix for the Hessian of the Lagrangian.",
            4: "The line search could not be terminated successfully.",
            5: "Length of a working array is too short. More detailed error information is obtained with IPRINT>0",
            6: "There are false dimensions, for example M>MMAX, N>=NMAX, or MNN2<>M+N+N+2.",
            7: "The search direction is close to zero, but the current iterate is still infeasible.",
            8: "The starting point violates a lower or upper bound.",
            9: "Wrong input parameter, i.e., MODE, LDL decomposition in D and C (in case of MODE=1), IPRINT, IOUT",
            10: "Internal inconsistency of the quadratic subproblem, division by zero.",
            11: "More than MAXFUN successive non-evaluable function calls.",
            100: (
                "The solution of the quadratic programming subproblem has been"
                + " terminated with an error message and IFAIL is set to IFQL+100,"
                + " where IFQL denotes the index of an inconsistent constraint."
            ),
        }
        return informs

    @staticmethod
    def _getDefaultOptions():
        defOpts = {
            "accuracy": [float, 1e-6],
            "accuracyQP": [float, 1e-14],
            "stepMin": [float, 1e-6],
            "maxFun": [int, 20],
            "maxIt": [int, 500],
            "maxNM": [int, 1],
            "rho": [float, 1.0],
            "iPrint": [int, [2, 0, 1, 3, 4]],
            "mode": [int, 0],
            "iOut": [int, 6],
            "lMerit": [bool, True],
            "lQl": [bool, False],
            "iFile": [str, "NLPQLP.out"],
        }
        return defOpts

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
            Specify method to compute sensitivities. To
            explicitly use pyOptSparse gradient class to do the
            derivatives with finite differences use 'FD'. 'sens'
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
            optimization.  The optimization problem used to generate
            the history file specified in 'hotStart' must be
            **IDENTICAL** to the currently supplied 'optProb'. By
            identical we mean, **EVERY SINGLE PARAMETER MUST BE
            IDENTICAL**. As soon as he requested evaluation point
            from NLPQL does not match the history, function and
            gradient evaluations revert back to normal evaluations.

        storeSens : bool
            Flag specifying if sensitivities are to be stored in hist.
            This is necessary for hot-starting only.
        """
        self.startTime = time.time()
        self.callCounter = 0
        self.storeSens = storeSens

        if len(optProb.constraints) == 0:
            self.unconstrained = True
            optProb.dummyConstraint = True

        # Save the optimization problem and finalize constraint
        # Jacobian, in general can only do on root proc
        self.optProb = optProb
        self.optProb.finalize()
        # Set history/hotstart/coldstart
        self._setHistory(storeHistory, hotStart)
        self._setInitialCacheValues()
        self._setSens(sens, sensStep, sensMode)
        blx, bux, xs = self._assembleContinuousVariables()
        xs = np.maximum(xs, blx)
        xs = np.minimum(xs, bux)
        nvar = len(xs)
        f = self._assembleObjective()

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
            # NLPQL - Objective/Constraint Values Function (Real Valued)
            # =================================================================
            def nlfunc(m, me, mmax, n, f, g, x, active, fail):
                fobj, fcon, fail = self._masterFunc(x, ["fobj", "fcon"])
                f = fobj
                g[0:m] = -fcon
                return f, g, fail

            # =================================================================
            # NLPQL - Objective/Constraint Gradients Function
            # =================================================================
            def nlgrad(m, me, mmax, n, f, g, df, dg, x, active, wa):
                gobj, gcon, fail = self._masterFunc(x, ["gobj", "gcon"])
                df[0:n] = gobj.copy()
                dg[0:m, 0:n] = -gcon.copy()
                return df, dg

            # setup argument list values

            num_procs = 1  # We only allow a single "processor" ie we are
            # actually running NLPQL (no P)

            # Set som basic sizes
            m = m
            me = meq
            mmax = max(1, m)

            n = nvar
            nmax = max(2, n + 2)
            mnn2 = m + n + n + 2

            # xs, ff, and gg have to have an extra dimension
            # associated with them for the NP. We will do this
            # correctly even though num_procs is hard-coded to 1.
            xs = np.array(xs).T
            f = np.array(f)
            g = np.zeros((mmax, num_procs))

            df = np.zeros(nmax)
            dg = np.zeros((mmax, nmax))
            u = np.zeros(mnn2)
            c = np.zeros((nmax, nmax))
            d = np.zeros(nmax)
            go = self.getOption
            if go("iPrint") < 0 or go("iPrint") > 4:
                raise ValueError("Incorrect iPrint option. Must be >=0 and <= 4")

            if not (go("mode") >= 0 and go("mode") <= 18):
                raise ValueError("Incorrect mode option. Must be >= 0 and <= 18.")

            if os.path.isfile(go("iFile")):
                os.remove(go("iFile"))
            ifail = np.array(0, dtype=int)
            # Run NLPQL
            t0 = time.time()
            # fmt: off
            nlpqlp.wrapper(num_procs, m, me, mmax, n, nmax, mnn2, xs, f, g, df, dg, u,
                           blx, bux, c, d, go('accuracy'), go('accuracyQP'),
                           go('stepMin'), go('maxFun'), go('maxIt'), go('maxNM'),
                           go('rho'), go('mode'), go('iPrint'), go('iOut'),
                           go('iFile'), ifail, go('lMerit'), go('lQl'),
                           nlfunc, nlgrad)
            # fmt: on
            optTime = time.time() - t0

            # Broadcast a -1 to indcate NLPQL has finished
            self.optProb.comm.bcast(-1, root=0)

            if self.storeHistory:
                self.metadata["endTime"] = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                self.metadata["optTime"] = optTime
                self.hist.writeData("metadata", self.metadata)
                self.hist.close()

            # Store Results
            inform = ifail.item()
            sol_inform = {}
            sol_inform["value"] = inform
            sol_inform["text"] = self.informs[inform]

            # Create the optimization solution
            sol = self._createSolution(optTime, sol_inform, f, xs)

        else:  # We are not on the root process so go into waiting loop:
            self._waitLoop()
            sol = None

        # Communication solution and return
        sol = self._communicateSolution(sol)

        return sol
