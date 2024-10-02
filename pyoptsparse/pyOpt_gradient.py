# Standard Python modules
from typing import Tuple, Union

# External modules
import numpy as np
from numpy import ndarray

# Local modules
from .pyOpt_MPI import MPI
from .pyOpt_optimization import Optimization
from .pyOpt_types import Dict1DType, Dict2DType


class Gradient:
    def __init__(self, optProb: Optimization, sensType: str, sensStep: float = None, sensMode: str = "", comm=None):
        """
        Gradient class for automatically computing gradients with finite
        difference or complex step.

        Parameters
        ----------
        optProb : Optimization instance
            This is the complete description of the optimization problem.

        sensType : str
            - ``FD`` for forward difference
            - ``CD`` for central difference
            - ``FDR`` for forward difference with relative step size
            - ``CDR`` for central difference with relative step size
            - ``CS`` for complex step

        sensStep : float
            Step size to use for differencing

        sensMode : str
            Flag to compute gradients in parallel.
        """
        self.optProb = optProb
        self.sensType = sensType
        self.sensStep: Union[float, complex]
        if sensStep is None:
            if self.sensType in ["fd", "fdr"]:
                self.sensStep = 1e-6
            elif self.sensType in ["cd", "cdr"]:
                self.sensStep = 1e-4
            else:
                self.sensStep = 1e-40j
        else:
            self.sensStep = sensStep
        self.sensMode = sensMode
        self.comm = comm

        # Now we can compute which dvs each process will need to
        # compute:
        ndvs = self.optProb.ndvs
        if self.sensMode == "pgc" and self.comm:
            self.mydvs = list(range(self.comm.rank, ndvs, self.comm.size))
        else:
            self.mydvs = list(range(ndvs))

    def _eval_func(self, x: ndarray) -> Tuple[ndarray, ndarray, bool]:
        """internal method to call function and extract obj, con"""

        xCall = self.optProb.processXtoDict(x)
        # Call objective
        [funcs, fail] = self.optProb.objFun(xCall)

        # Process constraint in case they are in dict form
        self.optProb.evaluateLinearConstraints(x, funcs)
        fobj = self.optProb.processObjtoVec(funcs, scaled=False)

        if self.sensType == "cs":
            fcon = self.optProb.processContoVec(funcs, scaled=False, dtype="D", natural=True)
        else:
            fcon = self.optProb.processContoVec(funcs, scaled=False, natural=True)

        return fobj, fcon, fail

    def __call__(self, x: Dict1DType, funcs: Dict1DType) -> Tuple[Dict2DType, bool]:
        """
        We need to make this object "look" the same as a user supplied
        function handle. That way, the optimizers need not care how
        the gradients are **actually** calculated.

        Parameters
        ----------
        x : dict
            Optimization variables from optimizer

        funcs : dict
            Dictionary of all the function values

        Returns
        -------
        funcsSens : dict
            Dictionary of sensitivities

        fail : bool
            Flag for failure. It currently always returns False
        """

        # Since this is *very* dumb loop over all the design
        # variables, it is easier to just loop over the x values as an
        # array. Furthermore, since the user **should** have
        # reasonably well scaled variables, the fixed step size should
        # have more meaning.

        # Generate final array sizes for the objective and constraint
        # gradients
        ndvs = self.optProb.ndvs
        ncon = self.optProb.nCon
        gobj = np.zeros(ndvs, "d")
        gcon = np.zeros((ncon, ndvs), "d")

        if self.sensMode == "pgc":
            funcsBase = self.comm.bcast(funcs)
        else:
            funcsBase = funcs

        # We DO NOT want the constraints scaled here....the constraint
        # scaling will be taken into account when the derivatives are
        # processed as per normal.
        xBase = self.optProb.processXtoVec(x)
        self.optProb.evaluateLinearConstraints(xBase, funcsBase)
        fconBase = self.optProb.processContoVec(funcsBase, scaled=False, natural=True)
        fobjBase = self.optProb.processObjtoVec(funcsBase, scaled=False)

        # Convert to complex if necessary:
        if self.sensType == "cs":
            xBase = xBase.astype("D")

        masterFail = False

        for i in self.mydvs:
            xph = xBase.copy()
            if self.sensType in ["fdr", "cdr"]:
                sensStep = max(abs(self.sensStep * xph[i]), self.sensStep)
            else:
                sensStep = self.sensStep
            xph[i] += sensStep

            fobj_ph, fcon_ph, fail = self._eval_func(xph)
            if fail:
                masterFail = True

            # forward difference
            if self.sensType in ["fd", "fdr"]:
                gobj[i] = (fobj_ph - fobjBase) / sensStep
                gcon[:, i] = (fcon_ph - fconBase) / sensStep

            # central difference
            elif self.sensType in ["cd", "cdr"]:
                xmh = xph  # xph isn't used anymore so  just point to same location
                xmh[i] -= 2 * sensStep

                fobj_mh, fcon_mh, fail = self._eval_func(xmh)
                if fail:
                    masterFail = True

                gobj[i] = (fobj_ph - fobj_mh) / (2 * sensStep)
                gcon[:, i] = (fcon_ph - fcon_mh) / (2 * sensStep)

            # complex step
            else:
                gobj[i] = np.imag(fobj_ph) / np.imag(sensStep)
                gcon[:, i] = np.imag(fcon_ph) / np.imag(sensStep)

        if self.sensMode == "pgc":
            # We just mpi_reduce to the root with sum. This uses the
            # efficient numpy versions
            self.comm.Reduce(gobj.copy(), gobj, op=MPI.SUM, root=0)
            self.comm.Reduce(gcon.copy(), gcon, op=MPI.SUM, root=0)

        # Logically reduce (over the comm) if the fail if *ANY*
        # gradient calc failed:
        if self.comm is not None:
            masterFail = self.comm.allreduce(masterFail, op=MPI.LOR)

        # Finally, we have to convert everything **back** to a
        # dictionary so the rest of the code works:
        funcsSens: Dict2DType = {}
        for objKey in self.optProb.objectives:
            funcsSens[objKey] = {}
            for dvGroup in self.optProb.variables:
                ss = self.optProb.dvOffset[dvGroup]
                funcsSens[objKey][dvGroup] = gobj[ss[0] : ss[1]]

        for conKey in self.optProb.constraints:
            con = self.optProb.constraints[conKey]
            funcsSens[conKey] = {}
            for dvGroup in self.optProb.variables:
                ss = self.optProb.dvOffset[dvGroup]
                funcsSens[conKey][dvGroup] = gcon[con.rs : con.re, ss[0] : ss[1]]

        return funcsSens, masterFail
