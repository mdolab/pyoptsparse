# /bin/env python
"""
pySNOPT - A variation of the pySNOPT wrapper specificially designed to
work with sparse optimization problems.
"""
# =============================================================================
# SNOPT Library
# =============================================================================
try:
    from . import snopt
except ImportError:
    snopt = None
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

# # ===========================================================================
# # Extension modules
# # ===========================================================================
from ..pyOpt_optimizer import Optimizer
from ..pyOpt_error import Error
from ..pyOpt_utils import ICOL, IDATA, IROW, extractRows, mapToCSC, scaleRows

# =============================================================================
# SNOPT Optimizer Class
# =============================================================================
class SNOPT(Optimizer):
    """
    SNOPT Optimizer Class - Inherited from Optimizer Abstract Class
    """

    def __init__(self, raiseError=True, *args, **kwargs):
        """
        SNOPT Optimizer Class Initialization
        """

        name = "SNOPT"
        category = "Local Optimizer"
        self.defOpts = {
            # SNOPT Printing Options
            "Major print level": [int, 1],  # Majors Print (1 - line major iteration log)
            "Minor print level": [int, 1],  # Minors Print (1 - line minor iteration log)
            "Print file": [str, "SNOPT_print.out"],  # Print File Name (specified by subroutine snInit)
            "iPrint": [int, 18],  # Print File Output Unit (override internally in snopt?)
            "Summary file": [str, "SNOPT_summary.out"],  # Summary File Name (specified by subroutine snInit)
            "iSumm": [int, 19],  # Summary File Output Unit (override internally in snopt?)
            "Print frequency": [int, 100],  # Minors Log Frequency on Print File
            "Summary frequency": [int, 100],  # Minors Log Frequency on Summary File
            "Solution": [str, "Yes"],  # Print Solution on the Print File
            "Suppress options listing": [type(None), None],  # (options are normally listed)
            "System information": [str, "No"],  # Print System Information on the Print File
            # SNOPT Problem Specification Options
            "Problem Type": [str, "Minimize"],  # Or 'Maximize', or 'Feasible point'
            "Objective row": [int, 1],  # row number of objective in F(x) (has precedence over ObjRow (snOptA))
            "Infinite bound": [float, 1.0e20],  # Infinite Bound Value
            # SNOPT Convergence Tolerances Options
            "Major feasibility tolerance": [float, 1.0e-6],  # Target Nonlinear Constraint Violation
            "Major optimality tolerance": [float, 1.0e-6],  # Target Complementarity Gap
            "Minor feasibility tolerance": [float, 1.0e-6],  # For Satisfying the QP Bounds
            # SNOPT Derivative Checking Options
            "Verify level": [int, 0],  # Gradients Check Flag
            "Start objective check at column": [int, 1],  # Start the gradient verification at this column
            "Start constraint check at column": [int, 1],
            # SNOPT Scaling Options
            "Scale option": [int, 1],  # Scaling (1 - linear constraints and variables)
            "Scale tolerance": [float, 0.9],  # Scaling Tolerance
            "Scale Print": [type(None), None],  # Default: scales are not printed
            # SNOPT Other Tolerances Options
            "Crash tolerance": [float, 0.1],
            "Linesearch tolerance": [float, 0.9],  # smaller for more accurate search
            "Pivot tolerance": [float, 3.7e-11],  # epsilon^(2/3)
            # SNOPT QP subproblems Options
            "QPSolver": [str, "Cholesky"],  # Default: Cholesky
            "Crash option": [int, 3],  # (3 - first basis is essentially triangular)
            "Elastic mode": [str, "No"],  # (start with elastic mode until necessary)
            "Elastic weight": [float, 1.0e4],  # (used only during elastic mode)
            "Iterations limit": [int, 10000],  # (or 20*ncons if that is more)
            "Partial price": [int, 1],  # (10 for large LPs)
            "Start": [str, "Cold"],  # has precedence over argument start, ('Warm': alternative to a cold start)
            # SNOPT SQP method Options
            "Major iterations limit": [int, 1000],  # or ncons if that is more
            "Minor iterations limit": [int, 500],  # or 3*ncons if that is more
            "Major step limit": [float, 2.0],
            "Superbasics limit": [int, None],  # (n1 + 1, n1 = number of nonlinear variables)
            "Derivative level": [int, 3],  # (NOT ALLOWED IN snOptA)
            "Derivative option": [int, 1],  # (ONLY FOR snOptA)
            "Derivative linesearch": [type(None), None],
            "Nonderivative linesearch": [type(None), None],
            "Function precision": [float, 3.0e-13],  # epsilon^0.8 (almost full accuracy)
            "Difference interval": [float, 5.5e-7],  # Function precision^(1/2)
            "Central difference interval": [float, 6.7e-5],  # Function precision^(1/3)
            "New superbasics limit": [int, 99],  # controls early termination of QPs
            "Penalty parameter": [float, 0.0],  # initial penalty parameter
            "Proximal point method": [int, 1],  # (1 - satisfies linear constraints near x0)
            "Reduced Hessian dimension": [int, 2000],  # (or Superbasics limit if that is less)
            "Violation limit": [float, 10.0],  # (unscaled constraint violation limit)
            "Unbounded step size": [float, 1.0e18],
            "Unbounded objective": [float, 1.0e15],
            # SNOPT Hessian approximation Options
            "Hessian full memory": [type(None), None],  # default if n1 <= 75
            "Hessian limited memory": [type(None), None],  # default if n1 > 75
            "Hessian frequency": [int, 999999],  # for full Hessian (never reset)
            "Hessian updates": [int, 10],  # for limited memory Hessian
            "Hessian flush": [int, 999999],  # no flushing
            # SNOPT Frequencies Options
            "Check frequency": [int, 60],  # test row residuals ||Ax - sk||
            "Expand frequency": [int, 10000],  # for anti-cycling procedure
            "Factorization frequency": [int, 50],  # 100 for LPs
            "Save frequency": [int, 100],  # save basis map
            # SNOPT LUSOL Options
            "LU factor tolerance": [float, 3.99],  # for NP (100.0 for LP)
            "LU update tolerance": [float, 3.99],  # for NP ( 10.0 for LP)
            "LU singularity tolerance": [float, 3.2e-11],
            "LU partial pivoting": [type(None), None],  # default threshold pivoting strategy
            "LU rook pivoting": [type(None), None],  # threshold rook pivoting
            "LU complete pivoting": [type(None), None],  # threshold complete pivoting
            # SNOPT Basis files Options
            "Old basis file": [int, 0],  # input basis map
            "New basis file": [int, 0],  # output basis map
            "Backup basis file": [int, 0],  # output extra basis map
            "Insert file": [int, 0],  # input in industry format
            "Punch file": [int, 0],  # output Insert data
            "Load file": [int, 0],  # input names and values
            "Dump file": [int, 0],  # output Load data
            "Solution file": [int, 0],  # different from printed solution
            # SNOPT Partitions of cw, iw, rw Options
            "Total character workspace": [int, 500],  # lencw: 500
            "Total integer workspace": [int, None],  # leniw: 500 + 100 * (m+n)
            "Total real workspace": [int, None],  # lenrw: 500 + 200 * (m+n)
            "User character workspace": [int, 500],
            "User integer workspace": [int, 500],
            "User real workspace": [int, 500],
            # SNOPT Miscellaneous Options
            "Debug level": [int, 1],  # (0 - Normal, 1 - for developers)
            "Timing level": [int, 3],  # (3 - print cpu times)
            # pySNOPT Options
            "Save major iteration variables": [
                list,
                ["step", "merit", "feasibility", "optimality", "penalty"],
            ],  # 'Hessian', 'slack', 'lambda' and 'condZHZ' are also supported
        }
        self.informs = {
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

        if snopt is None:
            if raiseError:
                raise Error("There was an error importing the compiled snopt module")

        self.set_options = []
        Optimizer.__init__(self, name, category, self.defOpts, self.informs, *args, **kwargs)

        # SNOPT need Jacobians in csc format
        self.jacType = "csc"

        # SNOPT specific Jacobian map
        self._snopt_jac_map_csr_to_csc = None

        # Check if we have numpy version 1.13.1. This version broke the callback in snopt.
        if np.__version__ == "1.13.1":
            raise Error("SNOPT is not compatible with numpy 1.13.1. Please use a different numpy version")

    def __call__(
        self,
        optProb,
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
        self.optProb.finalizeDesignVariables()
        self.optProb.finalizeConstraints()

        self._setInitialCacheValues()
        self._setSens(sens, sensStep, sensMode)
        blx, bux, xs = self._assembleContinuousVariables()
        ff = self._assembleObjective()

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
                blc = [-1e20]
                buc = [1e20]

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
                    raise Error("Failed to properly open %s, ierror = %3d" % (PrintFile, ierror))

            iSumm = self.getOption("iSumm")
            SummFile = os.path.join(self.getOption("Summary file"))
            if iSumm != 0 and iSumm != 6:
                ierror = snopt.openunit(iSumm, SummFile, "replace", "sequential")
                if ierror != 0:
                    raise Error("Failed to properly open %s, ierror = %3d" % (SummFile, ierror))

            # Calculate the length of the work arrays
            # --------------------------------------
            nvar = self.optProb.ndvs
            lencw = 500
            leniw = 500 + 100 * (ncon + nvar)
            lenrw = 500 + 200 * (ncon + nvar)

            self.options["Total integer workspace"][1] = leniw
            self.options["Total real workspace"][1] = lenrw

            cw = np.empty((lencw, 8), "c")
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

            mincw, miniw, minrw, cw = snopt.snmemb(iExit, ncon, nvar, neA, neGcon, nnCon, nnJac, nnObj, cw, iw, rw)

            if (minrw > lenrw) or (miniw > leniw) or (mincw > lencw):
                if mincw > lencw:
                    lencw = mincw
                    cw = np.array((lencw, 8), "c")
                    cw[:] = " "
                if miniw > leniw:
                    leniw = miniw
                    iw = np.zeros(leniw, np.intc)
                if minrw > lenrw:
                    lenrw = minrw
                    rw = np.zeros(lenrw, np.float)

                snopt.sninit(iPrint, iSumm, cw, iw, rw)

                # snInit resets all the options to the defaults.
                # Set them again!
                self._set_snopt_options(iPrint, iSumm, cw, iw, rw)

            # Setup argument list values
            start = np.array(self.options["Start"][1])
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
            cu = np.array(["        "], "c")
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

            # Set history/hotstart
            self._setHistory(storeHistory, hotStart)

            # The snopt c interface
            timeA = time.time()
            # fmt: off
            snopt.snkerc(start, nnCon, nnObj, nnJac, iObj, ObjAdd, ProbNm,
                         self._userfg_wrap, snopt.snlog, snopt.snlog2, snopt.sqlog, self._snstop,
                         Acol, indA, locA, bl, bu, Names, hs, xs, pi, rc, inform,
                         mincw, miniw, minrw, nS, ninf, sinf, ff, cu, iu, ru, cw, iw, rw)
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
                snopt.closeunit(self.options["iPrint"][1])
            if iSumm != 0 and iSumm != 6:
                snopt.closeunit(self.options["iSumm"][1])

            # Store Results
            inform = np.asscalar(inform)
            sol_inform = {}
            sol_inform["value"] = inform
            sol_inform["text"] = self.informs[inform]

            # Create the optimization solution
            sol = self._createSolution(optTime, sol_inform, ff, xs[:nvar], multipliers=pi)

        else:  # We are not on the root process so go into waiting loop:
            self._waitLoop()
            sol = None

        # Communication solution and return
        sol = self._communicateSolution(sol)

        return sol

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

    def _set_snopt_options(self, iPrint, iSumm, cw, iw, rw):
        """
        Set all the options into SNOPT that have been assigned
        by the user
        """

        # Set Options from the local options dictionary
        # ---------------------------------------------
        inform = np.array([-1], np.intc)
        for item in self.set_options:
            name = item[0]
            value = item[1]

            if name == "iPrint" or name == "iSumm":
                continue

            if isinstance(value, str):
                if name == "Start":
                    if value == "Cold":
                        snopt.snset("Cold start", iPrint, iSumm, inform, cw, iw, rw)
                    elif value == "Warm":
                        snopt.snset("Warm start", iPrint, iSumm, inform, cw, iw, rw)
                elif name == "Problem Type":
                    if value == "Minimize":
                        snopt.snset("Minimize", iPrint, iSumm, inform, cw, iw, rw)
                    elif value == "Maximize":
                        snopt.snset("Maximize", iPrint, iSumm, inform, cw, iw, rw)
                    elif value == "Feasible point":
                        snopt.snset("Feasible point", iPrint, iSumm, inform, cw, iw, rw)
                elif name == "Print file":
                    snopt.snset(name + " " + "%d" % iPrint, iPrint, iSumm, inform, cw, iw, rw)
                elif name == "Summary file":
                    snopt.snset(name + " " + "%d" % iSumm, iPrint, iSumm, inform, cw, iw, rw)
                else:
                    snopt.snset(name + " " + value, iPrint, iSumm, inform, cw, iw, rw)
            elif isinstance(value, float):
                snopt.snsetr(name, value, iPrint, iSumm, inform, cw, iw, rw)
            elif isinstance(value, int):
                snopt.snseti(name, value, iPrint, iSumm, inform, cw, iw, rw)
            elif isinstance(value, type(None)):
                snopt.snset(name, iPrint, iSumm, inform, cw, iw, rw)

        return

    def _on_setOption(self, name, value):
        """
        Set Optimizer Option Value (Optimizer Specific Routine)

        Documentation last updated:  May. 07, 2008 - Ruben E. Perez
        """

        self.set_options.append([name, value])

    def _on_getOption(self, name):
        """
        Get Optimizer Option Value (Optimizer Specific Routine)

        Documentation last updated:  May. 07, 2008 - Ruben E. Perez
        """

        pass

    def _on_getInform(self, infocode):
        """
        Get Optimizer Result Information (Optimizer Specific Routine)

        Keyword arguments:
        -----------------
        id -> STRING: Option Name

        Documentation last updated:  May. 07, 2008 - Ruben E. Perez
        """

        #
        mjr_code = (infocode[0] / 10) * 10
        # mnr_code = infocode[0] - 10*mjr_code
        try:
            inform_text = self.informs[mjr_code]
        except KeyError:
            inform_text = "Unknown Exit Status"
        # end try

        return inform_text

    def _on_flushFiles(self):
        """
        Flush the Output Files (Optimizer Specific Routine)

        Documentation last updated:  August. 09, 2009 - Ruben E. Perez
        """

        #
        iPrint = self.options["iPrint"][1]
        iSumm = self.options["iSumm"][1]
        if iPrint != 0:
            snopt.pyflush(iPrint)

        if iSumm != 0:
            snopt.pyflush(iSumm)
