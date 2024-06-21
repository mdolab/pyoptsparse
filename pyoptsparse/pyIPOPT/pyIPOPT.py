"""
pyIPOPT - A python wrapper to the core IPOPT compiled module.
"""

# Standard Python modules
import copy
import datetime
import os
import time

# External modules
import numpy as np

# Local modules
from ..pyOpt_optimizer import Optimizer
from ..pyOpt_utils import (
    ICOL,
    INFINITY,
    IROW,
    convertToCOO,
    extractRows,
    scaleRows,
    try_import_compiled_module_from_path,
)

# import the compiled module
THIS_DIR = os.path.dirname(os.path.abspath(__file__))
pyipoptcore = try_import_compiled_module_from_path("pyipoptcore", THIS_DIR)


class IPOPT(Optimizer):
    """
    IPOPT Optimizer Class - Inherited from Optimizer Abstract Class
    """

    def __init__(self, raiseError=True, options={}):
        """
        IPOPT Optimizer Class Initialization
        """

        name = "IPOPT"
        category = "Local Optimizer"
        defOpts = self._getDefaultOptions()
        informs = self._getInforms()

        if isinstance(pyipoptcore, str) and raiseError:
            raise ImportError(pyipoptcore)

        super().__init__(
            name,
            category,
            defaultOptions=defOpts,
            informs=informs,
            options=options,
            checkDefaultOptions=False,
        )

        # IPOPT needs Jacobians in coo format
        self.jacType = "coo"

    @staticmethod
    def _getInforms():
        informs = {
            0: "Solve Succeeded",
            1: "Solved To Acceptable Level",
            2: "Infeasible Problem Detected",
            3: "Search Direction Becomes Too Small",
            4: "Diverging Iterates",
            5: "User Requested Stop",
            6: "Feasible Point Found",
            -1: "Maximum Iterations Exceeded",
            -2: "Restoration Failed",
            -3: "Error In Step Computation",
            -4: "Maximum CpuTime Exceeded",
            -10: "Not Enough Degrees Of Freedom",
            -11: "Invalid Problem Definition",
            -12: "Invalid Option",
            -13: "Invalid Number Detected",
            -100: "Unrecoverable Exception",
            -101: "NonIpopt Exception Thrown",
            -102: "Insufficient Memory",
            -199: "Internal Error",
        }
        return informs

    @staticmethod
    def _getDefaultOptions():
        defOpts = {
            "print_level": [int, 0],
            "file_print_level": [int, 5],
            "sb": [str, "yes"],
            "print_user_options": [str, "yes"],
            "output_file": [str, "IPOPT.out"],
            "linear_solver": [str, "mumps"],
        }
        return defOpts

    def __call__(
        self,
        optProb,
        sens=None,
        sensStep=None,
        sensMode=None,
        storeHistory=None,
        hotStart=None,
        storeSens=True,
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
            Specifiy method to compute sensitivities.  To explictly
            use pyOptSparse gradient class to do the derivatives with
            finite differenes use \'FD\'. \'sens\' may also be \'CS\'
            which will cause pyOptSpare to compute the derivatives
            using the complex step method. Finally, \'sens\' may be a
            python function handle which is expected to compute the
            sensitivities directly. For expensive function evaluations
            and/or problems with large numbers of design variables
            this is the preferred method.

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
            IDENTICAL**. As soon as he requested evaluation point does
            not match the history, function and gradient evaluations
            revert back to normal evaluations.

        storeSens : bool
            Flag sepcifying if sensitivities are to be stored in hist.
            This is necessay for hot-starting only.
        """
        self.startTime = time.time()
        self.callCounter = 0
        self.storeSens = storeSens

        self.userRequestedTermination = False

        if len(optProb.constraints) == 0:
            # If the user *actually* has an unconstrained problem,
            # IPOPT sort of chokes with that....it has to have at
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
        blx, bux, xs = self._assembleContinuousVariables()
        self._setSens(sens, sensStep, sensMode)

        # Determine the sparsity structure of the full Jacobian
        # -----------------------------------------------------

        # Gather dummy data and process Jacobian:
        gcon = {}
        for iCon in self.optProb.constraints:
            gcon[iCon] = self.optProb.constraints[iCon].jac

        jac = self.optProb.processConstraintJacobian(gcon)

        if self.optProb.nCon > 0:
            # We need to reorder this full Jacobian...so get ordering:
            indices, blc, buc, fact = self.optProb.getOrdering(["ne", "ni", "le", "li"], oneSided=False)
            self.optProb.jacIndices = indices
            self.optProb.fact = fact
            self.optProb.offset = np.zeros(len(indices))
            ncon = len(indices)
            jac = extractRows(jac, indices)  # Does reordering
            scaleRows(jac, fact)  # Perform logical scaling
        else:
            blc = np.array(-INFINITY)
            buc = np.array(INFINITY)
            ncon = 1

        jac = convertToCOO(jac)  # Conver to coo format for IPOPT

        # We make a split here: If the rank is zero we setup the
        # problem and run IPOPT, otherwise we go to the waiting loop:

        if self.optProb.comm.rank == 0:
            # Now what we need for IPOPT is precisely the .row and
            # .col attributes of the fullJacobian array
            matStruct = (
                jac["coo"][IROW].copy().astype("int_"),
                jac["coo"][ICOL].copy().astype("int_"),
            )

            # Define the 4 call back functions that ipopt needs:
            def eval_f(x, user_data=None):
                fobj, fail = self._masterFunc(x, ["fobj"])
                if fail == 1:
                    fobj = np.array(np.NaN)
                elif fail == 2:
                    self.userRequestedTermination = True
                return fobj

            def eval_g(x, user_data=None):
                fcon, fail = self._masterFunc(x, ["fcon"])
                if fail == 1:
                    fcon = np.array(np.NaN)
                elif fail == 2:
                    self.userRequestedTermination = True
                return fcon.copy()

            def eval_grad_f(x, user_data=None):
                gobj, fail = self._masterFunc(x, ["gobj"])
                if fail == 1:
                    gobj = np.array(np.NaN)
                elif fail == 2:
                    self.userRequestedTermination = True
                return gobj.copy()

            def eval_jac_g(x, flag, user_data=None):
                if flag:
                    return copy.deepcopy(matStruct)
                else:
                    gcon, fail = self._masterFunc(x, ["gcon"])
                    if fail == 1:
                        gcon = np.array(np.NaN)
                    elif fail == 2:
                        self.userRequestedTermination = True
                    return gcon.copy()

            # Define intermediate callback. If this method returns false,
            # Ipopt will terminate with the User_Requested_Stop status.
            def eval_intermediate_callback(*args, **kwargs):
                if self.userRequestedTermination is True:
                    return False
                else:
                    return True

            timeA = time.time()
            nnzj = len(matStruct[0])
            nnzh = 0

            nlp = pyipoptcore.create(
                len(xs),
                blx,
                bux,
                ncon,
                blc,
                buc,
                nnzj,
                nnzh,
                eval_f,
                eval_grad_f,
                eval_g,
                eval_jac_g,
            )

            nlp.set_intermediate_callback(eval_intermediate_callback)

            self._set_ipopt_options(nlp)
            x, zl, zu, constraint_multipliers, obj, status = nlp.solve(xs)
            nlp.close()
            optTime = time.time() - timeA

            if self.storeHistory:
                self.metadata["endTime"] = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                self.metadata["optTime"] = optTime
                self.hist.writeData("metadata", self.metadata)
                self.hist.close()

            # Store Results
            sol_inform = {}
            sol_inform["value"] = status
            sol_inform["text"] = self.informs[status]

            # Create the optimization solution
            sol = self._createSolution(optTime, sol_inform, obj, x)

            # Indicate solution finished
            self.optProb.comm.bcast(-1, root=0)
        else:
            self._waitLoop()
            sol = None

        # Communication solution and return
        sol = self._communicateSolution(sol)

        return sol

    def _set_ipopt_options(self, nlp):
        """
        set all of the the options in self.options in the ipopt instance nlp
        """
        # Set Options from the local options dictionary
        # ---------------------------------------------

        for name, value in self.options.items():
            if isinstance(value, str):
                nlp.str_option(name, value)
            elif isinstance(value, float):
                nlp.num_option(name, value)
            elif isinstance(value, int):
                nlp.int_option(name, value)
            else:
                print("invalid option type", type(value))
