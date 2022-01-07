"""
pySNOPT - A variation of the pySNOPT wrapper specificially designed to
work with sparse optimization problems.
"""
# Compiled module
try:
    from . import snopt  # isort: skip
except ImportError:
    snopt = None
# Standard Python modules
import datetime
import os
import re
import time
from typing import Any, Dict, Optional, Tuple

# External modules
from baseclasses.utils import CaseInsensitiveSet
import numpy as np
from numpy import ndarray

# Local modules
from ..pyOpt_error import Error
from ..pyOpt_optimization import Optimization
from ..pyOpt_optimizer import Optimizer
from ..pyOpt_utils import ICOL, IDATA, INFINITY, IROW, extractRows, mapToCSC, scaleRows


class SNOPT(Optimizer):
    """
    SNOPT Optimizer Class - Inherited from Optimizer Abstract Class
    """

    def __init__(self, raiseError=True, options: Dict = {}):
        """
        SNOPT Optimizer Class Initialization
        """

        name = "SNOPT"
        category = "Local Optimizer"
        defOpts = self._getDefaultOptions()
        # these are SNOPT-related options that do not get set via snset
        self.specialOptions = CaseInsensitiveSet(
            {
                "iPrint",
                "iSumm",
                "Start",
            }
        )
        # this is purely within pySNOPT, nothing to do with SNOPT itself
        self.pythonOptions = CaseInsensitiveSet({"Save major iteration variables"})

        informs = self._getInforms()

        if snopt is None:
            if raiseError:
                raise Error("There was an error importing the compiled snopt module")
            else:
                version = None
        else:
            # extract SNOPT version
            version_str = snopt.sntitle().decode("utf-8")
            # The version_str is going to look like
            # S N O P T  7.7.5    (Oct 2020)
            # we search between "S N O P T" and "("
            res = re.search(r"S N O P T(.*)\(", version_str)
            if res is not None:
                version = res.group(1).strip()
            else:
                version = None

        super().__init__(
            name,
            category,
            defaultOptions=defOpts,
            informs=informs,
            options=options,
            checkDefaultOptions=False,
            version=version,
        )

        # SNOPT need Jacobians in csc format
        self.jacType = "csc"

        # SNOPT specific Jacobian map
        self._snopt_jac_map_csr_to_csc: Optional[Tuple[ndarray, ndarray, ndarray]] = None

    @staticmethod
    def _getDefaultOptions() -> Dict[str, Any]:
        defOpts = {
            "iPrint": [int, 18],
            "iSumm": [int, 19],
            "Print file": [str, "SNOPT_print.out"],
            "Summary file": [str, "SNOPT_summary.out"],
            "Problem Type": [str, ["Minimize", "Maximize", "Feasible point"]],
            "Start": [str, ["Cold", "Warm"]],
            "Derivative level": [int, 3],
            "Proximal iterations limit": [int, 10000],
            "Total character workspace": [int, None],
            "Total integer workspace": [int, None],
            "Total real workspace": [int, None],
            "Save major iteration variables": [list, ["step", "merit", "feasibility", "optimality", "penalty"]],
        }
        return defOpts

    @staticmethod
    def _getInforms() -> Dict[int, str]:
        informs = {
            0: "finished successfully",
            1: "optimality conditions satisfied",
            2: "feasible point found",
            3: "requested accuracy could not be achieved",
            4: "weak QP minimizer",
            10: "the problem appears to be infeasible",
            11: "infeasible linear constraints",
            12: "infeasible linear equalities",
            13: "nonlinear infeasibilities minimized",
            14: "infeasibilities minimized",
            15: "infeasible linear constraints in QP subproblem",
            20: "the problem appears to be unbounded",
            21: "unbounded objective",
            22: "constraint violation limit reached",
            30: "resource limit error",
            31: "iteration limit reached",
            32: "major iteration limit reached",
            33: "the superbasics limit is too small",
            40: "terminated after numerical difficulties",
            41: "current point cannot be improved",
            42: "singular basis",
            43: "cannot satisfy the general constraints",
            44: "ill-conditioned null-space basis",
            50: "error in the user-supplied functions",
            51: "incorrect objective  derivatives",
            52: "incorrect constraint derivatives",
            53: "the QP Hessian is indefinite",
            54: "incorrect second derivatives",
            55: "incorrect derivatives",
            56: "irregular or badly scaled problem functions",
            60: "undefined user-supplied functions",
            61: "undefined function at the first feasible point",
            62: "undefined function at the initial point",
            63: "unable to proceed into undefined region",
            70: "user requested termination",
            71: "terminated during function evaluation",
            72: "terminated during constraint evaluation",
            73: "terminated during objective evaluation",
            74: "terminated from monitor routine",
            80: "insufficient storage allocated",
            81: "work arrays must have at least 500 elements",
            82: "not enough character storage",
            83: "not enough integer storage",
            84: "not enough real storage",
            90: "input arguments out of range",
            91: "invalid input argument",
            92: "basis file dimensions do not match this problem",
            93: "the QP Hessian is indefinite",
            100: "finished successfully",
            101: "SPECS file read",
            102: "Jacobian structure estimated",
            103: "MPS file read",
            104: "memory requirements estimated",
            105: "user-supplied derivatives appear to be correct",
            106: "no derivatives were checked",
            107: "some SPECS keywords were not recognized",
            110: "errors while processing MPS data",
            111: "no MPS file specified",
            112: "problem-size estimates too small",
            113: "fatal error in the MPS file",
            120: "errors while estimating Jacobian structure",
            121: "cannot find Jacobian structure at given point",
            130: "fatal errors while reading the SP",
            131: "no SPECS file (iSpecs le 0 or iSpecs gt 99)",
            132: "End-of-file while looking for a BEGIN",
            133: "End-of-file while reading SPECS file",
            134: "ENDRUN found before any valid SPECS",
            140: "system error",
            141: "wrong no of basic variables",
            142: "error in basis package",
        }
        return informs

    def __call__(
        self,
        optProb: Optimization,
        sens=None,
        sensStep=None,
        sensMode=None,
        storeHistory=None,
        hotStart=None,
        storeSens=True,
        timeLimit=None,
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
            Specifiy method to compute sensitivities. The default is
            None which will use SNOPT's own finite differences which
            are vastly superiour to the pyOptSparse implementation. To
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
            from SNOPT does not match the history, function and
            gradient evaluations revert back to normal evaluations.

        storeSens : bool
            Flag sepcifying if sensitivities are to be stored in hist.
            This is necessay for hot-starting only.

        timeLimit : float
            Specify the maximum amount of time for optimizer to run.
            Must be in seconds. This can be useful on queue systems when
            you want an optimization to cleanly finish before the
            job runs out of time.
        """

        self.callCounter = 0
        self.storeSens = storeSens

        # Store the starting time if the keyword timeLimit is given:
        self.timeLimit = timeLimit
        self.startTime = time.time()

        if len(optProb.constraints) == 0:
            # If the user *actually* has an unconstrained problem,
            # snopt sort of chokes with that....it has to have at
            # least one constraint. So we will add one
            # automatically here:
            self.unconstrained = True
            optProb.dummyConstraint = True

        self.optProb = optProb
        self.optProb.finalize()
        # Set history/hotstart
        self._setHistory(storeHistory, hotStart)
        self._setInitialCacheValues()
        self._setSens(sens, sensStep, sensMode)
        blx, bux, xs = self._assembleContinuousVariables()

        oneSided = False
        # Set the number of nonlinear constraints snopt *thinks* we have:
        if self.unconstrained:
            nnCon = 1
        else:
            indices, tmp1, tmp2, fact = self.optProb.getOrdering(["ne", "ni"], oneSided=oneSided)
            nnCon = len(indices)
            self.optProb.jacIndices = indices
            self.optProb.fact = fact
            self.optProb.offset = np.zeros_like(fact)

            # Again, make SNOPT think we have a nonlinear constraint when all
            # our constraints are linear
            if nnCon == 0:
                nnCon = 1
                self.optProb.jacIndices = [0]
                self.optProb.fact = np.array([1.0])
                self.optProb.offset = np.zeros_like(self.optProb.fact)
        sol = None
        # We make a split here: If the rank is zero we setup the
        # problem and run SNOPT, otherwise we go to the waiting loop:
        if self.optProb.comm.rank == 0:

            # Determine the sparsity structure of the full Jacobian
            # -----------------------------------------------------

            # Gather dummy data and process Jacobian:
            gcon = {}
            for iCon in self.optProb.constraints:
                gcon[iCon] = self.optProb.constraints[iCon].jac

            jac = self.optProb.processConstraintJacobian(gcon)

            if self.optProb.nCon > 0:
                # We need to reorder this full Jacobian...so get ordering:
                indices, blc, buc, fact = self.optProb.getOrdering(["ne", "ni", "le", "li"], oneSided=oneSided)
                jac = extractRows(jac, indices)  # Does reordering
                scaleRows(jac, fact)  # Perform logical scaling
            else:
                blc = [-INFINITY]
                buc = [INFINITY]

            if self._snopt_jac_map_csr_to_csc is None:
                self._snopt_jac_map_csr_to_csc = mapToCSC(jac)

            # # CSC data is the csr data with the csc_indexing applied
            Acol = jac["csr"][IDATA][self._snopt_jac_map_csr_to_csc[IDATA]]
            # # CSC Row indices are just the row indices information from the map
            indA = self._snopt_jac_map_csr_to_csc[IROW] + 1
            # # CSC Column pointers are the column information from the map
            locA = self._snopt_jac_map_csr_to_csc[ICOL] + 1

            if self.optProb.nCon == 0:
                ncon = 1
            else:
                ncon = len(indices)

            # Initialize the Print and Summary files
            # --------------------------------------
            iPrint = self.getOption("iPrint")
            PrintFile = os.path.join(self.getOption("Print file"))
            if iPrint != 0 and iPrint != 6:
                ierror = snopt.openunit(iPrint, PrintFile, "replace", "sequential")
                if ierror != 0:
                    raise Error(f"Failed to properly open {PrintFile}, ierror = {ierror:3}")

            iSumm = self.getOption("iSumm")
            SummFile = os.path.join(self.getOption("Summary file"))
            if iSumm != 0 and iSumm != 6:
                ierror = snopt.openunit(iSumm, SummFile, "replace", "sequential")
                if ierror != 0:
                    raise Error(f"Failed to properly open {SummFile}, ierror = {ierror:3}")

            # Calculate the length of the work arrays
            # ---------------------------------------
            nvar = self.optProb.ndvs
            lencw = self.getOption("Total character workspace")
            leniw = self.getOption("Total integer workspace")
            lenrw = self.getOption("Total real workspace")

            # Set flags to avoid overwriting user-specified lengths
            checkLencw = lencw is None
            checkLeniw = leniw is None
            checkLenrw = lenrw is None

            # Set defaults
            minWorkArrayLength = 500
            if lencw is None:
                lencw = minWorkArrayLength
                self.setOption("Total character workspace", lencw)
            if leniw is None:
                leniw = minWorkArrayLength + 100 * (ncon + nvar)
                self.setOption("Total integer workspace", leniw)
            if lenrw is None:
                lenrw = minWorkArrayLength + 200 * (ncon + nvar)
                self.setOption("Total real workspace", lenrw)

            cw = np.empty((lencw, 8), dtype="|S1")
            iw = np.zeros(leniw, np.intc)
            rw = np.zeros(lenrw, np.float)
            snopt.sninit(iPrint, iSumm, cw, iw, rw)

            # Memory allocation
            nnObj = nvar
            nnJac = nvar
            iObj = np.array(0, np.intc)
            neA = len(indA)
            neGcon = neA  # The nonlinear Jacobian and A are the same
            iExit = 0

            # Set the options into the SNOPT instance
            self._set_snopt_options(iPrint, iSumm, cw, iw, rw)

            # Estimate workspace storage requirement
            mincw, miniw, minrw, cw = snopt.snmemb(iExit, ncon, nvar, neA, neGcon, nnCon, nnJac, nnObj, cw, iw, rw)

            # This flag is set to True if any of the lengths are overwritten
            lengthsChanged = False

            # Overwrite lengths if the defaults are too small
            if checkLencw and mincw > lencw:
                lencw = mincw
                self.setOption("Total character workspace", lencw)
                cw = np.empty((lencw, 8), dtype="|S1")
                cw[:] = " "
                lengthsChanged = True
            if checkLeniw and miniw > leniw:
                leniw = miniw
                self.setOption("Total integer workspace", leniw)
                iw = np.zeros(leniw, np.intc)
                lengthsChanged = True
            if checkLenrw and minrw > lenrw:
                lenrw = minrw
                self.setOption("Total real workspace", lenrw)
                rw = np.zeros(lenrw, np.float)
                lengthsChanged = True

            # Initialize SNOPT again if any of the lengths were overwritten
            if lengthsChanged:
                snopt.sninit(iPrint, iSumm, cw, iw, rw)

                # snInit resets all the options to the defaults.
                # Set them again!
                self._set_snopt_options(iPrint, iSumm, cw, iw, rw)

            # Setup argument list values
            start = np.array(self.getOption("Start"))
            ObjAdd = np.array([0.0], np.float)
            ProbNm = np.array(self.optProb.name, "c")
            cdummy = -1111111  # this is a magic variable defined in SNOPT for undefined strings
            cw[51, :] = cdummy  # we set these to cdummy so that a placeholder is used in printout
            cw[52, :] = cdummy
            cw[53, :] = cdummy
            cw[54, :] = cdummy
            xs = np.concatenate((xs, np.zeros(ncon, np.float)))
            bl = np.concatenate((blx, blc))
            bu = np.concatenate((bux, buc))
            leniu = 2
            lenru = 3
            cu = np.empty((1, 8), dtype="|S1")
            iu = np.zeros(leniu, np.intc)
            ru = np.zeros(lenru, np.float)
            hs = np.zeros(nvar + ncon, np.intc)

            Names = np.array(["        "], "c")
            pi = np.zeros(ncon, np.float)
            rc = np.zeros(nvar + ncon, np.float)
            inform = np.array([-1], np.intc)
            mincw = np.array([0], np.intc)
            miniw = np.array([0], np.intc)
            minrw = np.array([0], np.intc)
            nS = np.array([0], np.intc)
            ninf = np.array([0], np.intc)
            sinf = np.array([0.0], np.float)

            # The snopt c interface
            timeA = time.time()
            # fmt: off
            obj = snopt.snkerc(start, nnCon, nnObj, nnJac, iObj, ObjAdd, ProbNm,
                               self._userfg_wrap, snopt.snlog, snopt.snlog2, snopt.sqlog, self._snstop,
                               Acol, indA, locA, bl, bu, Names, hs, xs, pi, rc, inform,
                               mincw, miniw, minrw, nS, ninf, sinf, cu, iu, ru, cw, iw, rw)
            # fmt: on
            optTime = time.time() - timeA

            # Indicate solution finished
            self.optProb.comm.bcast(-1, root=0)

            if self.storeHistory:
                # Record the full state of variables, xs and hs such
                # that we could perform a warm start.
                self.hist.writeData("xs", xs)
                self.hist.writeData("hs", hs)
                self.metadata["endTime"] = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                self.metadata["optTime"] = optTime
                self.hist.writeData("metadata", self.metadata)
                self.hist.close()

            if iPrint != 0 and iPrint != 6:
                snopt.closeunit(self.getOption("iPrint"))
            if iSumm != 0 and iSumm != 6:
                snopt.closeunit(self.getOption("iSumm"))

            # Store Results
            inform = inform.item()
            sol_inform = {}
            sol_inform["value"] = inform
            sol_inform["text"] = self.informs[inform]

            # Create the optimization solution
            sol = self._createSolution(optTime, sol_inform, obj, xs[:nvar], multipliers=pi)

        else:  # We are not on the root process so go into waiting loop:
            self._waitLoop()

        # Communication solution and return
        commSol = self._communicateSolution(sol)

        return commSol

    def _userfg_wrap(self, mode, nnJac, x, fobj, gobj, fcon, gcon, nState, cu, iu, ru):
        """
        The snopt user function. This is what is actually called from snopt.

        Essentially nothing is done in this function, but this funcion
        has to precisely match the signature from fortran so must look
        EXACTLY like this.

        All we do here is call the generic masterFunc in the baseclass
        which will take care of everything else.
        """
        # nState >=2 means this is the final call which is redundant
        # here we just return without doing anything since we don't
        # need to do any cleanup or anything
        if nState >= 2:
            return

        fail = 0
        if mode == 0 or mode == 2:
            fobj, fcon, fail = self._masterFunc(x, ["fobj", "fcon"])
        if fail == 0:
            if mode == 1:
                if self.getOption("Derivative level") != 0:
                    gobj, gcon, fail = self._masterFunc(x, ["gobj", "gcon"])
            if mode == 2:
                if self.getOption("Derivative level") != 0:
                    gobj, gcon, fail2 = self._masterFunc(x, ["gobj", "gcon"])
                    fail = max(fail, fail2)

        if fail == 1:
            mode = -1
        elif fail == 2:
            mode = -2

        # Flush the files to the buffer for all the people who like to
        # monitor the residual
        snopt.pyflush(self.getOption("iPrint"))
        snopt.pyflush(self.getOption("iSumm"))

        # Check if we've exceeded the timeLimit
        if self.timeLimit is not None:
            if time.time() - self.startTime > self.timeLimit:
                mode = -2  # User requested termination

        return mode, fobj, gobj, fcon, gcon

    def _getHessian(self, iw, rw):
        """
        This function retrieves the approximate Hessian from the SNOPT workspace arrays
        Call it for example from the _snstop routine or after SNOPT has finished, where iw and rw arrays are available
        Currently only full memory Hessian mode is implemented, do not use this for limited-memory case.

        The FM Hessian in SNOPT is stored with its Cholesky factor
        which has been flattened to 1D
        """
        lvlHes = iw[72 - 1]  # 0,1,2 => LM, FM, Exact Hessian
        if lvlHes != 1:
            print("pyOptSparse Error! Limited-memory Hessian not supported for history file!")
            return None
        lU = iw[391 - 1] - 1  # U(lenU), BFGS Hessian H = U'U
        lenU = iw[392 - 1]
        Uvec = rw[lU : lU + lenU]
        nnH = iw[24 - 1]
        Umat = np.zeros((nnH, nnH))
        Umat[np.triu_indices(nnH)] = Uvec
        H = np.matmul(Umat.T, Umat)
        return H

    def _getPenaltyParam(self, iw, rw):
        """
        Retrieves the full penalty parameter vector from the work arrays.
        """
        nnCon = iw[23 - 1]
        lxPen = iw[304 - 1] - 1
        xPen = rw[lxPen : lxPen + nnCon]
        return xPen

    # fmt: off
    def _snstop(self, ktcond, mjrprtlvl, minimize, n, nncon, nnobj, ns, itn, nmajor, nminor, nswap, condzhz,
                iobj, scaleobj, objadd, fobj, fmerit, penparm, step, primalinf, dualinf, maxvi, maxvirel, hs,
                locj, indj, jcol, scales, bl, bu, fx, fcon, gcon, gobj, ycon, pi, rc, rg, x, cu, iu, ru, cw, iw, rw):
    # fmt: on # noqa: E115
        """
        This routine is called every major iteration in SNOPT, after solving QP but before line search
        Currently we use it just to determine the correct major iteration counting,
        and save some parameters in history if needed

        returning with iabort != 0 will terminate SNOPT immediately
        """
        iterDict = {
            "isMajor": True,
            "nMajor": nmajor,
            "nMinor": nminor,
        }
        for saveVar in self.getOption("Save major iteration variables"):
            if saveVar == "merit":
                iterDict[saveVar] = fmerit
            elif saveVar == "feasibility":
                iterDict[saveVar] = primalinf
            elif saveVar == "optimality":
                iterDict[saveVar] = dualinf
            elif saveVar == "penalty":
                penParam = self._getPenaltyParam(iw, rw)
                iterDict[saveVar] = penParam
            elif saveVar == "Hessian":
                H = self._getHessian(iw, rw)
                iterDict[saveVar] = H
            elif saveVar == "step":
                iterDict[saveVar] = step
            elif saveVar == "condZHZ":
                iterDict[saveVar] = condzhz
            elif saveVar == "slack":
                iterDict[saveVar] = x[n:]
            elif saveVar == "lambda":
                iterDict[saveVar] = pi
        if self.storeHistory:
            currX = x[:n]  # only the first n component is x, the rest are the slacks
            if nmajor == 0:
                callCounter = 0
            else:
                xuser_vec = self.optProb._mapXtoUser(currX)
                callCounter = self.hist._searchCallCounter(xuser_vec)
            if callCounter is not None:
                self.hist.write(callCounter, iterDict)
        iabort = 0
        return iabort

    def _set_snopt_options(self, iPrint: int, iSumm: int, cw: ndarray, iw: ndarray, rw: ndarray):
        """
        Set all the options into SNOPT that have been assigned
        by the user
        """
        # Set Options from the local options dictionary
        # ---------------------------------------------
        inform = np.array([-1], np.intc)
        for name, value in self.options.items():
            # these do not get set using snset
            if name in self.specialOptions or name in self.pythonOptions:
                continue

            if isinstance(value, str):
                if name == "Problem Type":
                    snopt.snset(value, iPrint, iSumm, inform, cw, iw, rw)
                elif name == "Print file":
                    snopt.snset(name + " " + f"{iPrint}", iPrint, iSumm, inform, cw, iw, rw)
                elif name == "Summary file":
                    snopt.snset(name + " " + f"{iSumm}", iPrint, iSumm, inform, cw, iw, rw)
                else:
                    snopt.snset(name + " " + value, iPrint, iSumm, inform, cw, iw, rw)
            elif isinstance(value, float):
                snopt.snsetr(name, value, iPrint, iSumm, inform, cw, iw, rw)
            elif isinstance(value, int):
                snopt.snseti(name, value, iPrint, iSumm, inform, cw, iw, rw)
            elif isinstance(value, type(None)):
                snopt.snset(name, iPrint, iSumm, inform, cw, iw, rw)

    def _on_flushFiles(self):
        """
        Flush the Output Files (Optimizer Specific Routine)
        """

        #
        iPrint = self.getOption("iPrint")
        iSumm = self.getOption("iSumm")
        if iPrint != 0:
            snopt.pyflush(iPrint)

        if iSumm != 0:
            snopt.pyflush(iSumm)
