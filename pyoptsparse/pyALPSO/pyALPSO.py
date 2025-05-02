"""
pyALPSO - A pyOptSparse interface to ALPSO
work with sparse optimization problems.
"""

# Standard Python modules
import datetime
import time

# External modules
import numpy as np

# Local modules
from . import alpso
from ..pyOpt_optimizer import Optimizer

# isort: off


class ALPSO(Optimizer):
    """
    ALPSO Optimizer Class - Inherited from Optimizer Abstract Class

    *Keyword arguments:**

    - pll_type -> STR: ALPSO Parallel Implementation (None, SPM- Static, DPM- Dynamic, POA-Parallel Analysis), *Default* = None
    """

    def __init__(self, options={}):
        self.alpso = alpso

        category = "Global Optimizer"
        defOpts = self._getDefaultOptions()
        informs = self._getInforms()
        super().__init__("ALPSO", category, defaultOptions=defOpts, informs=informs, options=options)

    @staticmethod
    def _getInforms():
        informs = {}
        return informs

    @staticmethod
    def _getDefaultOptions():
        defOpts = {
            "SwarmSize": [int, 40],
            "maxOuterIter": [int, 200],
            "maxInnerIter": [int, 6],
            "minInnerIter": [int, 6],
            "dynInnerIter": [int, 0],
            "stopCriteria": [int, 1],
            "stopIters": [int, 5],
            "etol": [float, 1e-3],
            "itol": [float, 1e-3],
            # 'ltol':[float, 1e-2],
            "rtol": [float, 1e-2],
            "atol": [float, 1e-2],
            "dtol": [float, 1e-1],
            "printOuterIters": [int, 0],
            "printInnerIters": [int, 0],
            "rinit": [float, 1.0],
            "xinit": [int, 0],
            "vinit": [float, 1.0],
            "vmax": [float, 2.0],
            "c1": [float, 2.0],
            "c2": [float, 1.0],
            "w1": [float, 0.99],
            "w2": [float, 0.55],
            "ns": [int, 15],
            "nf": [int, 5],
            "dt": [float, 1.0],
            "vcrazy": [float, 1e-4],
            "fileout": [int, 1],
            "filename": [str, "ALPSO.out"],
            "seed": [int, 0],
            "HoodSize": [int, 40],
            "HoodModel": [str, "gbest"],
            "HoodSelf": [int, 1],
            "Scaling": [int, 1],
            "parallelType": [str, [None, "EXT"]],
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
            optimization.  The optimization problem used to generate
            the history file specified in 'hotStart' must be
            **IDENTICAL** to the currently supplied 'optProb'. By
            identical we mean, **EVERY SINGLE PARAMETER MUST BE
            IDENTICAL**. As soon as he requested evaluation point
            from ALPSO does not match the history and function
            evaluations revert back to normal evaluations.

        Notes
        -----
        The kwargs are there such that the sens= argument can be
        supplied (but ignored here in alpso)
        """
        self.startTime = time.time()
        # ======================================================================
        # ALPSO - Objective/Constraint Values Function
        # ======================================================================

        def objconfunc(x):
            fobj, fcon, fail = self._masterFunc(x, ["fobj", "fcon"])
            return fobj, fcon

        # Save the optimization problem and finalize constraint
        # Jacobian, in general can only do on root proc
        self.optProb = optProb
        self.optProb.finalize()
        # Set history/hotstart/coldstart
        self._setHistory(storeHistory, hotStart)
        self._setInitialCacheValues()

        if len(optProb.constraints) == 0:
            self.unconstrained = True

        xl, xu, xs = self._assembleContinuousVariables()
        xs = np.maximum(xs, xl)
        xs = np.minimum(xs, xu)
        n = len(xs)
        types = [0] * len(xs)
        oneSided = True
        if self.unconstrained:
            m = 0
            me = 0
        else:
            indices, blc, buc, fact = self.optProb.getOrdering(
                ["ne", "le", "ni", "li"], oneSided=oneSided, noEquality=False
            )
            m = len(indices)

            self.optProb.jacIndices = indices
            self.optProb.fact = fact
            self.optProb.offset = buc
            indices, __, __, __ = self.optProb.getOrdering(["ne", "le"], oneSided=oneSided, noEquality=False)
            me = len(indices)

        if self.optProb.comm.rank == 0:
            # Setup argument list values
            opt = self.getOption

            dyniI = self.getOption("dynInnerIter")
            if dyniI == 0:
                self.setOption("minInnerIter", opt("maxInnerIter"))

            if opt("stopCriteria") not in [0, 1]:
                raise ValueError("Incorrect Stopping Criteria Setting")

            if opt("fileout") not in [0, 1, 2, 3]:
                raise ValueError("Incorrect fileout Setting")

            # Run ALPSO
            t0 = time.time()
            # fmt: off
            opt_x, opt_f, opt_g, opt_lambda, nfevals, rseed = self.alpso.alpso(
                n, m, me, types, xs, xl, xu, opt('SwarmSize'), opt('HoodSize'),
                opt('HoodModel'), opt('maxOuterIter'), opt('maxInnerIter'),
                opt('minInnerIter'), opt('stopCriteria'), opt('stopIters'),
                opt('etol'), opt('itol'), opt('rtol'), opt('atol'), opt('dtol'),
                opt('printOuterIters'), opt('printInnerIters'), opt('rinit'),
                opt('vinit'), opt('vmax'), opt('c1'), opt('c2'), opt('w1'),
                opt('w2'), opt('ns'), opt('nf'), opt('vcrazy'), opt('fileout'),
                opt('filename'), None, None, opt('seed'),
                opt('Scaling'), opt('HoodSelf'), objconfunc)
            # fmt: on
            optTime = time.time() - t0

            if self.storeHistory:
                self.metadata["endTime"] = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                self.metadata["optTime"] = optTime
                self.hist.writeData("metadata", self.metadata)
                self.hist.close()

            # Broadcast a -1 to indcate NSGA2 has finished
            self.optProb.comm.bcast(-1, root=0)

            # Store Results
            sol_inform = {"value": "", "text": ""}

            # Create the optimization solution
            sol = self._createSolution(optTime, sol_inform, opt_f, opt_x)
        else:  # We are not on the root process so go into waiting loop:
            self._waitLoop()
            sol = None

        # Communication solution and return
        sol = self._communicateSolution(sol)

        return sol

    def _on_setOption(self, name, value):
        if name == "parallelType":
            if isinstance(value, str) and value.upper() == "EXT":
                try:
                    from . import alpso_ext

                    self.alpso = alpso_ext
                except ImportError:
                    raise ImportError("pyALPSO: ALPSO EXT shared library failed to import.")

    def _communicateSolution(self, sol):
        if sol is not None:
            sol.userObjCalls = self.optProb.comm.allreduce(sol.userObjCalls)
            sol.comm = None
        sol = self.optProb.comm.bcast(sol)
        sol.objFun = self.optProb.objFun
        sol.comm = self.optProb.comm

        return sol
