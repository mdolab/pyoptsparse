"""
pyCONMIN - A variation of the pyCONMIN wrapper specificially designed to
work with sparse optimization problems.
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
conmin = try_import_compiled_module_from_path("conmin", THIS_DIR, raise_warning=True)


class CONMIN(Optimizer):
    """
    CONMIN Optimizer Class - Inherited from Optimizer Abstract Class
    """

    def __init__(self, raiseError=True, options={}):
        name = "CONMIN"
        category = "Local Optimizer"
        defOpts = self._getDefaultOptions()
        informs = self._getInforms()
        if isinstance(conmin, str) and raiseError:
            raise ImportError(conmin)

        self.set_options = []
        super().__init__(name, category, defaultOptions=defOpts, informs=informs, options=options)

        # CONMIN needs Jacobians in dense format
        self.jacType = "dense2d"

    @staticmethod
    def _getInforms():
        informs = {}
        return informs

    @staticmethod
    def _getDefaultOptions():
        defOpts = {
            "ITMAX": [int, int(1e4)],
            "DELFUN": [float, 1e-6],
            "DABFUN": [float, 1e-6],
            "ITRM": [int, 5],
            "NFEASCT": [int, 20],
            "IPRINT": [int, 4],
            "IOUT": [int, 6],
            "IFILE": [str, "CONMIN.out"],
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
            from CONMIN does not match the history, function and
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
        # Set history/hotstart/coldstart
        self._setHistory(storeHistory, hotStart)
        self._setInitialCacheValues()
        self._setSens(sens, sensStep, sensMode)
        blx, bux, xs = self._assembleContinuousVariables()
        xs = np.maximum(xs, blx)
        xs = np.minimum(xs, bux)
        ff = self._assembleObjective()

        oneSided = True
        noEquality = True
        if self.unconstrained:
            m = 0
        else:
            indices, blc, buc, fact = self.optProb.getOrdering(
                ["ne", "le", "ni", "li"], oneSided=oneSided, noEquality=noEquality
            )
            m = len(indices)

            self.optProb.jacIndices = indices
            self.optProb.fact = fact
            self.optProb.offset = buc

        if self.optProb.comm.rank == 0:
            # =================================================================
            # CONMIN - Objective/Constraint Values Function
            # =================================================================
            def cnmnfun(n1, n2, x, f, g):
                fobj, fcon, fail = self._masterFunc(x[0:ndv], ["fobj", "fcon"])
                f = fobj
                g[0:ncn] = fcon

                return f, g

            # =================================================================
            # CONMIN - Objective/Constraint Gradients Function
            # =================================================================
            def cnmngrad(n1, n2, x, f, g, ct, df, a, ic, nac):
                gobj, gcon, fail = self._masterFunc(x[0:ndv], ["gobj", "gcon"])
                df[0:ndv] = gobj.copy()

                # Only assign the gradients for constraints that are
                # actually active:
                nac = 0
                for j in range(ncn):
                    if g[j] >= ct:
                        a[0:ndv, nac] = gcon[j, :]
                        ic[nac] = j + 1
                        nac += 1
                return df, a, ic, nac

            # Setup argument list values
            ndv = len(xs)
            ncn = m
            nn1 = ndv + 2
            nn2 = ncn + 2 * ndv
            nn3 = max(nn2, ndv)
            nn4 = max(nn2, ndv)
            nn5 = 2 * nn4

            if ncn > 0:
                gg = np.zeros(ncn, float)
            else:
                gg = np.array([0], float)

            if self.getOption("IPRINT") >= 0 and self.getOption("IPRINT") <= 4:
                iprint = self.getOption("IPRINT")
            else:
                raise ValueError("IPRINT option must be >= 0 and <= 4")

            iout = self.getOption("IOUT")
            ifile = self.getOption("IFILE")

            # Check if file exists and remove if necessary
            if iprint > 0:
                if os.path.isfile(ifile):
                    os.remove(ifile)

            itmax = self.getOption("ITMAX")
            delfun = self.getOption("DELFUN")

            # finit, ginit = cnmnfun([],[],xx,ff,gg)
            dabfun = self.getOption("DABFUN")

            itrm = self.getOption("ITRM")
            nfeasct = self.getOption("ITRM")
            nfdg = 1  # User will supply all gradients

            # Counters for functions and gradients
            nfun = 0
            ngrd = 0

            # Run CONMIN
            t0 = time.time()
            # fmt: off
            conmin.conmin(ndv, ncn, xs, blx, bux, ff, gg,
                          nn1, nn2, nn3, nn4, nn5,
                          iprint, iout, ifile, itmax, delfun, dabfun, itrm,
                          nfeasct, nfdg, nfun, ngrd, cnmnfun, cnmngrad)
            # fmt: on
            optTime = time.time() - t0

            if self.storeHistory:
                self.metadata["endTime"] = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                self.metadata["optTime"] = optTime
                self.hist.writeData("metadata", self.metadata)
                self.hist.close()

            if iprint > 0:
                conmin.closeunit(self.getOption("IOUT"))

            # Broadcast a -1 to indcate SLSQP has finished
            self.optProb.comm.bcast(-1, root=0)

            # Store Results
            sol_inform = {"value": "", "text": ""}

            # Create the optimization solution
            sol = self._createSolution(optTime, sol_inform, ff, xs)

        else:  # We are not on the root process so go into waiting loop:
            self._waitLoop()
            sol = None

        # Communication solution and return
        sol = self._communicateSolution(sol)

        return sol
