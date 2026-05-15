"""
pyUno - A Python wrapper to the Uno optimizer via unopy.
"""

# Standard Python modules
import datetime
import importlib.metadata as ilmd
import pathlib
import time
from typing import TextIO

# External modules
import numpy as np
from packaging.version import Version

# Local modules
from ..pyOpt_optimizer import Optimizer
from ..pyOpt_solution import SolutionInform
from ..pyOpt_utils import ICOL, INFINITY, IROW, convertToCOO, extractRows, import_module, scaleRows

unopy = import_module("unopy")


_UNOPY_MIN_VERSION = "0.4.7"


class _UserTerminationError(Exception):
    """Raised an Uno callback function returns 2."""

    pass


class Uno(Optimizer):
    """
    Uno Optimizer Class - Inherited from Optimizer Abstract Class.

    Uses the ``unopy`` Python bindings to the Uno unified nonlinear optimizer.
    """

    def __init__(self, raiseError=True, options=None):
        """
        Uno Optimizer Class Initialization.

        Parameters
        ----------
        raiseError : bool
            If True, raises an ImportError when unopy is not available.
        options : dict
            Dictionary of options to pass to the optimizer.
        """
        if options is None:
            options = {}
        name = "Uno"
        category = "Local Optimizer"
        defOpts = self._getDefaultOptions()
        informs = self._getInforms()

        if isinstance(unopy, Exception):
            if raiseError:
                raise unopy
        else:
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
        self.pythonOptions = {"preset", "logger_stream", "save_major_iteration_variables"}

        # Save the result object from the optimize call separately from the
        # pyoptsparse Solution object, in case the user wants more detail.
        self.result = None

        # Flag to track user-requested termination
        self._userRequestedTermination = False

        # Reference to logger file, if opened by this optimizer
        self._logger_file_stream = None

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
            # Special options handled differently by Uno
            "preset": [str, "filtersqp"],
            "logger_stream": [(TextIO, str, pathlib.Path), None],
            "save_major_iteration_variables": [list, []],
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
            "protected_actual_reduction_macheps_coefficient": [float, 10.0],
            "hessian_model": [str, "LBFGS"],
            "quasi_newton_memory_size": [int, 20],
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

        self._logger_file_stream = None
        self._userRequestedTermination = False

        if len(optProb.constraints) == 0:
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
                    raise _UserTerminationError("User requested termination.")
                return fobj

            def _constraints(x, con_val):
                fcon, fail = self._masterFunc(x, ["fcon"])
                if fail == 1:
                    raise ValueError("Constraint evaluation failed.")
                elif fail == 2:
                    self._userRequestedTermination = True
                    raise _UserTerminationError("User requested termination.")
                con_val[:] = fcon

            def _objective_gradient(x, grad):
                gobj, fail = self._masterFunc(x, ["gobj"])
                if fail == 1:
                    raise ValueError("Objective gradient evaluation failed.")
                elif fail == 2:
                    self._userRequestedTermination = True
                    raise _UserTerminationError("User requested termination.")
                grad[:] = gobj

            def _jacobian(x, jac_val):
                gcon_vals, fail = self._masterFunc(x, ["gcon"])
                if fail == 1:
                    raise ValueError("Constraint gradient evaluation failed.")
                elif fail == 2:
                    self._userRequestedTermination = True
                    raise _UserTerminationError("User requested termination.")
                jac_val[:] = gcon_vals

            # Major-iteration callback: fired once per accepted iterate by Uno's outer loop,
            # not on line-search / trust-region trial points. Mirrors _snstop in pySNOPT.
            nMajorCounter = [0]

            def _notifyAcceptableIterate(
                primals,
                lower_bound_multipliers,
                upper_bound_multipliers,
                constraint_multipliers,
                objective_multiplier,
                primal_feasibility_residual,
                stationarity_residual,
                complementarity_residual,
            ):
                iterDict = {
                    "isMajor": True,
                    "nMajor": nMajorCounter[0],
                    "primal_feasibility_residual": primal_feasibility_residual,
                    "stationarity_residual": stationarity_residual,
                    "complementarity_residual": complementarity_residual,
                    "objective_multiplier": objective_multiplier,
                }
                nMajorCounter[0] += 1

                # Opt-in heavier data. Callback args are read-only views into Uno memory,
                # so we copy before storing.
                for saveVar in self.getOption("save_major_iteration_variables"):
                    if saveVar == "lower_bound_multipliers":
                        iterDict[saveVar] = np.array(lower_bound_multipliers, copy=True)
                    elif saveVar == "upper_bound_multipliers":
                        iterDict[saveVar] = np.array(upper_bound_multipliers, copy=True)
                    elif saveVar == "constraint_multipliers":
                        iterDict[saveVar] = np.array(constraint_multipliers, copy=True)
                    else:
                        raise ValueError(
                            f"Received unknown Uno save variable '{saveVar}'. "
                            "Please see 'save_major_iteration_variables' option in the "
                            "pyOptSparse documentation under 'Uno'."
                        )

                if self.storeHistory:
                    # Slice to len(xs): some presets (e.g. ipopt) append slack variables
                    # to primals, making it longer than the original design variable vector.
                    xuser_vec = self.optProb._mapXtoUser(np.array(primals[: len(xs)], copy=True))
                    # Like IPOPT, Uno calls objective and constraints separately, so we find two call counters and append iter_dict to both counters.
                    call_counter_1 = self.hist._searchCallCounter(xuser_vec)
                    if call_counter_1 is None:
                        call_counter_2 = None
                    else:
                        call_counter_2 = self.hist._searchCallCounter(xuser_vec, last=call_counter_1 - 1)
                    for callCounter in [call_counter_2, call_counter_1]:
                        if callCounter is not None:
                            self.hist.write(callCounter, iterDict)

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

            try:
                self._set_uno_options(solver)
                solver.set_notify_acceptable_iterate_callback(_notifyAcceptableIterate)
                self.result = solver.optimize(model)
            except _UserTerminationError:
                # User requested termination during optimization
                # Create a result object indicating user termination
                pass
            finally:
                if self._logger_file_stream:
                    self._logger_file_stream.close()

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
        solver.set_preset(self.getOption("preset"))

        logger_stream = self.getOption("logger_stream")
        if logger_stream is not None:
            if isinstance(logger_stream, (str, pathlib.Path)):
                self._logger_file_stream = open(logger_stream, "w")
                solver.set_logger_stream(self._logger_file_stream)
            else:
                solver.set_logger_stream(logger_stream)

        for name, value in self.options.items():
            # skip the options we handle externally
            if name in self.pythonOptions:
                continue
            solver.set_option(name, value)
