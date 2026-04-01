"""
pyUno - A Python wrapper to the Uno optimizer via unopy.
"""

# Standard Python modules
import datetime
import importlib.metadata as ilmd
import time

# External modules
import numpy as np
from packaging.version import Version

# Local modules
from ..pyOpt_optimizer import Optimizer
from ..pyOpt_solution import SolutionInform
from ..pyOpt_utils import ICOL, INFINITY, IROW, convertToCOO, extractRows, import_module, scaleRows

unopy = import_module("unopy")


_UNOPY_MIN_VERSION = "0.4.0"


class Uno(Optimizer):
    """
    Uno Optimizer Class - Inherited from Optimizer Abstract Class.

    Uses the ``unopy`` Python bindings to the Uno unified nonlinear optimizer.
    """

    def __init__(self, raiseError=True, options={}):
        """
        Uno Optimizer Class Initialization.

        Parameters
        ----------
        raiseError : bool
            If True, raises an ImportError when unopy is not available.
        options : dict
            Dictionary of options to pass to the optimizer.
        """
        name = "Uno"
        category = "Local Optimizer"
        defOpts = self._getDefaultOptions()
        informs = self._getInforms()
        if isinstance(unopy, Exception) and raiseError:
            raise unopy

        unopy_version = ilmd.version("unopy")

        if Version(unopy_version) < Version(_UNOPY_MIN_VERSION):
            raise RuntimeError(
                "The pyoptsparse Uno interface requires unopy "
                f"{_UNOPY_MIN_VERSION} or later, but {unopy_version} is installed"
            )

        super().__init__(
            name,
            category,
            defaultOptions=defOpts,
            informs=informs,
            options=options,
            checkDefaultOptions=False,
        )

        # Uno needs Jacobians in COO format
        self.jacType = "coo"

        # Options handled outside of unopy's set_option interface.
        # 'preset' is applied via solver.set_preset() rather than set_option().
        self.pythonOptions = {"preset"}

        # Save the result object from the optimize call separately from the
        # pyoptsparse Solution object, in case the user wants more detail.
        self.result = None

        # Flag to track user-requested termination
        self._userRequestedTermination = False

    @staticmethod
    def _getInforms():
        """
        Return a dictionary of inform codes and their descriptions.

        Returns
        -------
        dict
            Mapping from integer inform code to description string.
        """
        informs = {
            0: "Success",
            1: "Iteration Limit Exceeded",
            2: "Time Limit Exceeded",
            3: "Evaluation Error",
            4: "Algorithmic Error",
            5: "User Termination",
        }
        return informs

    @staticmethod
    def _getDefaultOptions():
        """
        Return a dictionary of default options.

        Returns
        -------
        dict
            Default options with type and default value pairs.
        """
        defOpts = {
            "preset": [str, "filtersqp"],
            # Termination Options
            "primal_tolerance": [float, 1e-8],
            "dual_tolerance": [float, 1e-8],
            "loose_primal_tolerance": [float, 1e-6],
            "loose_dual_tolerance": [float, 1e-6],
            "loose_tolerance_iteration_threshold": [int, 15],
            "max_iterations": [int, 2000],
            "time_limit": [float, float("inf")],
            "unbounded_objective_threshold": [float, -1e20],
            # Logging and Output Options
            "logger": [str, "INFO"],
            "print_solution": [bool, False],
            "print_subproblem": [bool, False],
            "progress_norm": [str, "L1"],
            "residual_norm": [str, "INF"],
            # Numerical Options
            "residual_scaling_threshold": [float, 100.0],
            "protect_actual_reduction_against_roundoff": [bool, False],
            "protected_actual_reduction_macheps_coefficient": [float, 10.],
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
        This is the main routine used to solve the optimization problem.

        Parameters
        ----------
        optProb : Optimization or Solution class instance
            This is the complete description of the optimization problem
            to be solved by the optimizer.

        sens : str or python Function.
            Specify method to compute sensitivities. To explicitly
            use pyOptSparse gradient class to do the derivatives with
            finite differences use 'FD'. 'sens' may also be 'CS'
            which will cause pyOptSparse to compute the derivatives
            using the complex step method. Finally, 'sens' may be a
            python function handle which is expected to compute the
            sensitivities directly. For expensive function evaluations
            and/or problems with large numbers of design variables
            this is the preferred method.

        sensStep : float
            Set the step size to use for design variables. Defaults to
            1e-6 when sens is 'FD' and 1e-40j when sens is 'CS'.

        sensMode : str
            Use 'pgc' for parallel gradient computations. Only
            available with mpi4py and each objective evaluation is
            otherwise serial.

        storeHistory : str
            File name of the history file into which the history of
            this optimization will be stored.

        hotStart : str
            File name of the history file to "replay" for the
            optimization. The optimization problem used to generate
            the history file specified in 'hotStart' must be
            **IDENTICAL** to the currently supplied 'optProb'. By
            identical we mean, **EVERY SINGLE PARAMETER MUST BE
            IDENTICAL**. As soon as the requested evaluation point does
            not match the history, function and gradient evaluations
            revert back to normal evaluations.

        storeSens : bool
            Flag specifying if sensitivities are to be stored in hist.
            This is necessary for hot-starting only.
        """
        self.startTime = time.time()
        self.callCounter = 0
        self.storeSens = storeSens

        if len(optProb.constraints) == 0:
            # Uno requires at least one constraint, so add a dummy one
            # for unconstrained problems.
            self.unconstrained = True
            optProb.dummyConstraint = True

        # Save the optimization problem and finalize constraint Jacobian,
        # in general can only do on root proc
        self.optProb = optProb
        self.optProb.finalize()
        # Set history/hotstart
        self._setHistory(storeHistory, hotStart)
        self._setInitialCacheValues()
        blx, bux, xs = self._assembleContinuousVariables()
        self._setSens(sens, sensStep, sensMode)

        # Determine the sparsity structure of the full Jacobian
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
            blc = np.atleast_1d(-INFINITY)
            buc = np.atleast_1d(INFINITY)
            ncon = 1

        jac = convertToCOO(jac)  # Convert to COO format for Uno

        # We make a split here: If the rank is zero we setup the
        # problem and run Uno, otherwise we go to the waiting loop:
        if self.optProb.comm.rank == 0:
            row_indices = jac["coo"][IROW].copy().astype("int_")
            col_indices = jac["coo"][ICOL].copy().astype("int_")
            nnz = len(row_indices)

            # Define the four callback functions that Uno needs.
            # unopy now provides x and output arrays as numpy arrays directly.
            def _objective(x):
                fobj, fail = self._masterFunc(x, ["fobj"])
                if fail == 1:
                    raise ValueError("Objective evaluation failed.")
                elif fail == 2:
                    self._userRequestedTermination = True
                    raise KeyboardInterrupt("User requested termination.")
                return fobj

            def _constraints(x, con_val):
                fcon, fail = self._masterFunc(x, ["fcon"])
                if fail == 1:
                    raise ValueError("Constraint evaluation failed.")
                elif fail == 2:
                    self._userRequestedTermination = True
                    raise KeyboardInterrupt("User requested termination.")
                con_val[:] = fcon

            def _objective_gradient(x, grad):
                gobj, fail = self._masterFunc(x, ["gobj"])
                if fail == 1:
                    raise ValueError("Objective gradient evaluation failed.")
                elif fail == 2:
                    self._userRequestedTermination = True
                    raise KeyboardInterrupt("User requested termination.")
                grad[:] = gobj

            def _jacobian(x, jac_val):
                gcon_vals, fail = self._masterFunc(x, ["gcon"])
                if fail == 1:
                    raise ValueError("Constraint gradient evaluation failed.")
                elif fail == 2:
                    self._userRequestedTermination = True
                    raise KeyboardInterrupt("User requested termination.")
                jac_val[:] = gcon_vals

            timeA = time.time()

            model = unopy.Model(
                unopy.PROBLEM_NONLINEAR,
                len(xs),
                blx,
                bux,
                unopy.ZERO_BASED_INDEXING,
            )

            model.set_objective(unopy.MINIMIZE, _objective, _objective_gradient)

            model.set_constraints(
                ncon,
                _constraints,
                blc,
                buc,
                nnz,
                row_indices,
                col_indices,
                _jacobian,
            )

            model.set_initial_primal_iterate(xs)

            solver = unopy.UnoSolver()
            self._set_uno_options(solver)

            try:
                self.result = solver.optimize(model)
            except KeyboardInterrupt:
                # User requested termination during optimization
                # Create a result object indicating user termination
                pass

            optTime = time.time() - timeA

            if self.storeHistory:
                self.metadata["endTime"] = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                self.metadata["optTime"] = optTime
                self.hist.writeData("metadata", self.metadata)
                self.hist.close()

            # Map Uno optimization status to pyoptsparse inform code
            if self._userRequestedTermination:
                # User requested termination
                inform_code = 5
            else:
                inform_code = self.result.optimization_status.value
            sol_inform = SolutionInform.from_informs(self.informs, inform_code)

            # Extract only the original variables (some presets may add additional variables)
            if self.result is not None:
                x_sol = self.result.primal_solution[: len(xs)]
                multipliers = self.result.constraint_dual_solution
                obj_val = self.result.solution_objective
            else:
                # No result available due to user termination
                x_sol = xs.copy()
                multipliers = None
                obj_val = np.nan

            # Create the optimization solution
            sol = self._createSolution(optTime, sol_inform, obj_val, x_sol, multipliers=multipliers)

            # Indicate solution finished
            self.optProb.comm.bcast(-1, root=0)
        else:
            self._waitLoop()
            sol = None

        # Communicate solution and return
        sol = self._communicateSolution(sol)

        return sol

    def _set_uno_options(self, solver):
        """
        Set all options in self.options on the Uno solver instance.

        Parameters
        ----------
        solver : unopy.UnoSolver
            The Uno solver instance to configure.
        """
        preset = self.getOption("preset")
        solver.set_preset(preset)

        for name, value in self.options.items():
            # skip preset
            if name in self.pythonOptions:
                continue
            solver.set_option(name, value)
