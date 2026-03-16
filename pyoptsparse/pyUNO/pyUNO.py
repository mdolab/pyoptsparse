"""
pyUNO - A Python wrapper to the UNO optimizer via unopy.
"""

# Standard Python modules
import datetime
import importlib.metadata as ilmd
import time
from packaging.version import Version

# External modules
import numpy as np

# Local modules
from ..pyOpt_optimizer import Optimizer
from ..pyOpt_solution import SolutionInform
from ..pyOpt_utils import ICOL, INFINITY, IROW, convertToCOO, extractRows, import_module, scaleRows

unopy = import_module('unopy')

# Mapping from unopy OptimizationStatus integer value to pyoptsparse inform codes.
# These values correspond to: SUCCESS=0, ITERATION_LIMIT=1, TIME_LIMIT=2,
# EVALUATION_ERROR=3, ALGORITHMIC_ERROR=4, UNBOUNDED=6.
_UNO_STATUS_TO_INFORM = {
    0: 0,    # SUCCESS
    1: -1,   # ITERATION_LIMIT
    2: -2,   # TIME_LIMIT
    3: -3,   # EVALUATION_ERROR
    4: -4,   # ALGORITHMIC_ERROR
    6: -5,   # UNBOUNDED
}


class UNO(Optimizer):
    """
    UNO Optimizer Class - Inherited from Optimizer Abstract Class.

    Uses the ``unopy`` Python bindings to the UNO unified nonlinear optimizer.
    UNO supports presets ``ipopt``, and ``filterslp``.
    """

    def __init__(self, raiseError=True, options={}):
        """
        UNO Optimizer Class Initialization.

        Parameters
        ----------
        raiseError : bool
            If True, raises an ImportError when unopy is not available.
        options : dict
            Dictionary of options to pass to the optimizer.
        """
        name = 'UNO'
        category = 'Local Optimizer'
        defOpts = self._getDefaultOptions()
        informs = self._getInforms()
        if isinstance(unopy, Exception) and raiseError:
            raise unopy
        
        if Version(ilmd.version('unopy')) < Version('0.3.0'):
            raise RuntimeError('The pyoptsparse UNO interface requires version 0.3.0 or later.')

        super().__init__(
            name,
            category,
            defaultOptions=defOpts,
            informs=informs,
            options=options,
            checkDefaultOptions=False,
        )

        # UNO needs Jacobians in COO format
        self.jacType = 'coo'

        # Options handled outside of unopy's set_option interface.
        # 'preset' is applied via solver.set_preset() rather than set_option().
        self.pythonOptions = {'preset'}

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
            0: 'Optimal',
            -1: 'Iteration Limit Exceeded',
            -2: 'Time Limit Exceeded',
            -3: 'Evaluation Error',
            -4: 'Algorithmic Error',
            -5: 'Unbounded',
            -99: 'Unknown Status',
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
            'preset': [str, 'filtersqp'],
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
            # UNO requires at least one constraint, so add a dummy one
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
            indices, blc, buc, fact = self.optProb.getOrdering(
                ['ne', 'ni', 'le', 'li'], oneSided=False
            )
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

        jac = convertToCOO(jac)  # Convert to COO format for UNO

        # We make a split here: If the rank is zero we setup the
        # problem and run UNO, otherwise we go to the waiting loop:
        if self.optProb.comm.rank == 0:
            row_indices = jac['coo'][IROW].copy().astype('int_')
            col_indices = jac['coo'][ICOL].copy().astype('int_')
            nnz = len(row_indices)

            # Define the four callback functions that UNO needs.
            # unopy now provides x and output arrays as numpy arrays directly.
            def _objective(x):
                fobj, fail = self._masterFunc(x, ['fobj'])
                if fail:
                    raise ValueError('Objective evaluation failed.')
                return fobj

            def _constraints(x, con_val):
                fcon, fail = self._masterFunc(x, ['fcon'])
                if fail:
                    raise ValueError('Constraint evaluation failed.')
                con_val[:] = fcon

            def _objective_gradient(x, grad):
                gobj, fail = self._masterFunc(x, ['gobj'])
                if fail:
                    raise ValueError('Objective gradient evaluation failed.')
                grad[:] = gobj

            def _jacobian(x, jac_val):
                gcon_vals, fail = self._masterFunc(x, ['gcon'])
                if fail:
                    raise ValueError('Constraint gradient evaluation failed.')
                jac_val[:] = gcon_vals

            timeA = time.time()

            model = unopy.Model(
                unopy.PROBLEM_NONLINEAR,
                len(xs),
                blx.tolist(),
                bux.tolist(),
                unopy.ZERO_BASED_INDEXING,
            )
            model.set_objective(unopy.MINIMIZE, _objective, _objective_gradient)
            model.set_constraints(
                ncon,
                _constraints,
                blc.tolist(),
                buc.tolist(),
                nnz,
                row_indices.tolist(),
                col_indices.tolist(),
                _jacobian,
            )
            model.set_initial_primal_iterate(xs.tolist())

            solver = unopy.UnoSolver()
            self._set_uno_options(solver)
            result = solver.optimize(model)

            optTime = time.time() - timeA

            if self.storeHistory:
                self.metadata['endTime'] = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
                self.metadata['optTime'] = optTime
                self.hist.writeData('metadata', self.metadata)
                self.hist.close()

            # Map UNO optimization status to pyoptsparse inform code
            status_val = result.optimization_status.value
            inform_code = _UNO_STATUS_TO_INFORM.get(status_val, -99)
            sol_inform = SolutionInform.from_informs(self.informs, inform_code)

            x_sol = np.array(list(result.primal_solution))
            multipliers = np.array(list(result.constraint_dual_solution))

            # Create the optimization solution
            sol = self._createSolution(
                optTime, sol_inform, result.solution_objective, x_sol, multipliers=multipliers
            )

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
        Set all options in self.options on the UNO solver instance.

        Parameters
        ----------
        solver : unopy.UnoSolver
            The UNO solver instance to configure.
        """
        solver.set_preset(self.getOption('preset'))
        for name, value in self.options.items():
            # skip pyUNO-specific options
            if name in self.pythonOptions:
                continue
            solver.set_option(name, str(value))
