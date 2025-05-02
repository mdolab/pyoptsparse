"""
pyNSGA2 - A variation of the pyNSGA2 wrapper specificially designed to
work with sparse optimization problems.
"""

# Standard Python modules
import os
import time

# External modules
import numpy as np

# Local modules
from ..pyOpt_optimizer import Optimizer
from ..pyOpt_utils import try_import_compiled_module_from_path

# import the compiled module
THIS_DIR = os.path.dirname(os.path.abspath(__file__))
nsga2 = try_import_compiled_module_from_path("nsga2", THIS_DIR, raise_warning=True)


class NSGA2(Optimizer):
    """
    NSGA2 Optimizer Class - Inherited from Optimizer Abstract Class
    """

    def __init__(self, raiseError=True, options={}):
        name = "NSGA-II"
        category = "Global Optimizer"
        defOpts = self._getDefaultOptions()
        informs = self._getInforms()
        super().__init__(name, category, defaultOptions=defOpts, informs=informs, options=options)

        if isinstance(nsga2, str) and raiseError:
            raise ImportError(nsga2)

        if self.getOption("PopSize") % 4 != 0:
            raise ValueError("Option 'PopSize' must be a multiple of 4")

    @staticmethod
    def _getInforms():
        informs = {}
        return informs

    @staticmethod
    def _getDefaultOptions():
        defOpts = {
            "PopSize": [int, 100],
            "maxGen": [int, 1000],
            "pCross_real": [float, 0.6],
            "pMut_real": [float, 0.2],
            "eta_c": [float, 10.0],
            "eta_m": [float, 20.0],
            "pCross_bin": [float, 0.0],
            "pMut_bin": [float, 0.0],
            "PrintOut": [int, 1],
            "seed": [int, 0],
            "xinit": [int, 0],
        }
        return defOpts

    def __call__(self, optProb, storeHistory=None, hotStart=None, **kwargs):
        """
        This is the main routine used to solve the optimization
        problem.

        Parameters
        ----------
        optProb : Optimization or Solution class instance
            This is the complete description of the optimization problem
            to be solved by the optimizer

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
            from NSGA2 does not match the history, function and
            gradient evaluations revert back to normal evaluations.

        Notes
        -----
        The kwargs are there such that the sens= argument can be
        supplied (but ignored here in nsga2)
        """

        # ======================================================================
        # NSGA-II - Objective/Constraint Values Function
        # ======================================================================
        def objconfunc(nreal, nobj, ncon, x, f, g):
            xx = np.array(x)
            fobj, fcon, fail = self._masterFunc(xx, ["fobj", "fcon"])
            fobj = np.atleast_1d(fobj)
            f[0:nobj] = fobj
            g[0:ncon] = -fcon[0:ncon]

            return f, g

        self.startTime = time.time()
        self.callCounter = 0

        if len(optProb.constraints) == 0:
            # If the user *actually* has an unconstrained problem,
            # slsqp sort of chokes with that....it has to have at
            # least one constraint. So we will add one
            # automatically here:
            self.unconstrained = True
            optProb.dummyConstraint = False

        # Save the optimization problem and finalize constraint
        # Jacobian, in general can only do on root proc
        self.optProb = optProb
        self.optProb.finalize()
        # Set history/hotstart
        self._setHistory(storeHistory, hotStart)
        self._setInitialCacheValues()

        blx, bux, xs = self._assembleContinuousVariables()
        xs = np.maximum(xs, blx)
        xs = np.minimum(xs, bux)
        n = len(xs)
        ff = self._assembleObjective()
        oneSided = True

        # Set the number of nonlinear constraints snopt *thinks* we have:
        if self.unconstrained:
            m = 0
        else:
            indices, blc, buc, fact = self.optProb.getOrdering(
                ["ne", "le", "ni", "li"], oneSided=oneSided, noEquality=True
            )
            m = len(indices)

            self.optProb.jacIndices = indices
            self.optProb.fact = fact
            self.optProb.offset = buc

        g = nsga2.new_doubleArray(m)
        len_ff = len(np.atleast_1d(ff))
        f = nsga2.new_doubleArray(len_ff)

        if self.optProb.comm.rank == 0:
            # Variables Handling
            n = len(xs)
            x = nsga2.new_doubleArray(n)
            xl = nsga2.new_doubleArray(n)
            xu = nsga2.new_doubleArray(n)
            for i in range(n):
                nsga2.doubleArray_setitem(x, i, xs[i])
                nsga2.doubleArray_setitem(xl, i, blx[i])
                nsga2.doubleArray_setitem(xu, i, bux[i])

            # Setup argument list values
            nfeval = 0
            opt = self.getOption

            if self.getOption("PrintOut") >= 0 and self.getOption("PrintOut") <= 2:
                printout = self.getOption("PrintOut")
            else:
                raise ValueError("Incorrect option PrintOut")

            seed = self.getOption("seed")
            if seed == 0:
                seed = time.time()

            # Run NSGA-II
            nsga2.set_pyfunc(objconfunc)
            t0 = time.time()
            # fmt: off
            nsga2.nsga2(n, m, len_ff, f, x, g, nfeval, xl, xu, opt('PopSize'), opt('maxGen'),
                        opt('pCross_real'), opt('pMut_real'), opt('eta_c'), opt('eta_m'),
                        opt('pCross_bin'), opt('pMut_bin'), printout, seed, opt('xinit'))
            # fmt: on
            optTime = time.time() - t0

            # Broadcast a -1 to indcate NSGA2 has finished
            self.optProb.comm.bcast(-1, root=0)

            # Store Results
            sol_inform = {"value": "", "text": ""}

            xstar = [0.0] * n
            for i in range(n):
                xstar[i] = nsga2.doubleArray_getitem(x, i)

            fStar = np.zeros(len_ff)
            if len_ff > 1:
                for i in range(len_ff):
                    fStar[i] = nsga2.doubleArray_getitem(f, i)
            else:
                fStar = nsga2.doubleArray_getitem(f, 0)

            # Create the optimization solution
            sol = self._createSolution(optTime, sol_inform, fStar, xstar)

        else:  # We are not on the root process so go into waiting loop:
            self._waitLoop()
            sol = None

        # Communication solution and return
        sol = self._communicateSolution(sol)

        return sol
