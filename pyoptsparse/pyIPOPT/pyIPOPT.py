# /bin/env python
"""
pyIPOPT - A python wrapper to the core IPOPT compiled module.
"""
# =============================================================================
# IPOPT Library
# =============================================================================

try:
    from . import pyipoptcore
except ImportError:
    pyipoptcore = None

# =============================================================================
# standard Python modules
# =============================================================================
import copy
import time
import datetime

# =============================================================================
# External Python modules
# =============================================================================
import numpy as np

# =============================================================================
# Extension modules
# =============================================================================
from ..pyOpt_optimizer import Optimizer
from ..pyOpt_error import Error
from ..pyOpt_utils import IROW, ICOL, convertToCOO, extractRows, scaleRows

# =============================================================================
# IPOPT Optimizer Class
# =============================================================================
class IPOPT(Optimizer):
    """
    IPOPT Optimizer Class - Inherited from Optimizer Abstract Class
    """

    def __init__(self, raiseError=True, *args, **kwargs):
        """
        IPOPT Optimizer Class Initialization
        """

        name = "IPOPT"
        category = "Local Optimizer"

        # These options are documented on the website:
        # http://www.coin-or.org/Ipopt/documentation/node39.html
        # accessed on March 26, 2014.

        self.defOpts = {
            "print_level": [int, 0],
            "output_file": [str, "IPOPT.out"],
            # Output verbosity level. '0-12'
            "file_print_level": [int, 5],
            "print_user_options": [str, "no"],
            "print_options_documentation": [str, "no"],
            "print_timing_statistics": [str, "no"],
            "option_file_name": [str, "IPOPT_options.opt"],
            "replace_bounds": [str, "no"],
            "skip_finalize_solution_call": [str, "no"],
            # Enables printing of additional info string at
            # end of iteration output.
            "print_info_string": [str, "no"],
            # Determines what value is printed in the "inf_pr"
            # output column. 'internal' or 'original'
            "inf_pr_output": [str, "original"],  #
            "print_frequency_iter": [int, 1],
            # Determines at which time frequency the
            # summarizing iteration output line should be
            # printed.
            "print_frequency_time": [float, 0.0],  #
            # Convergence.
            # 's_max' : 100.0, #
            "tol": [float, 1e-06],
            "max_iter": [int, 100],
            "max_cpu_time": [float, 1e06],
            "dual_inf_tol": [float, 1.0],
            "constr_viol_tol": [float, 0.0001],
            "compl_inf_tol": [float, 0.0001],
            "acceptable_tol": [float, 1e-06],
            "acceptable_iter": [int, 15],
            "acceptable_dual_inf_tol": [float, 1e10],
            "acceptable_constr_viol_tol": [float, 0.01],
            "acceptable_compl_inf_tol": [float, 0.01],
            "acceptable_obj_change_tol": [float, 1e20],
            "diverging_iterates_tol": [float, 1e20],
            "mu_target": [float, 0.0],
            # NLP Scaling.
            "nlp_scaling_method": [str, "gradient-based"],
            "obj_scaling_factor": [float, 1.0],
            "nlp_scaling_max_gradient": [float, 100.0],
            "nlp_scaling_obj_target_gradient": [float, 0.0],
            "nlp_scaling_constr_target_gradient": [float, 0.0],
            "nlp_scaling_min_value": [float, 1e-08],
            # NLP.
            "nlp_lower_bound_inf": [float, -1e19],
            "nlp_upper_bound_inf": [float, 1e19],
            "fixed_variable_treatment": [str, "make_parameter"],
            "dependency_detector": [str, "none"],
            "dependency_detection_with_rhs": [str, "no"],
            "num_linear_variables": [int, 0],
            "kappa_d": [float, 1e-05],
            "bound_relax_factor": [float, 1e-08],
            "honor_original_bounds": [str, "yes"],
            "check_derivatives_for_naninf": [str, "no"],
            "jac_c_constant": [str, "no"],
            "jac_d_constant": [str, "no"],
            "hessian_constant": [str, "no"],
            # Initialization.
            "bound_push": [float, 0.01],
            "bound_frac": [float, 0.01],
            "slack_bound_push": [float, 0.01],
            "slack_bound_frac": [float, 0.01],
            "constr_mult_init_max": [float, 1000.0],
            "bound_mult_init_val": [float, 1.0],
            "bound_mult_init_method": [str, "constant"],
            "least_square_init_primal": [str, "no"],
            "least_square_init_duals": [str, "no"],
            # Barrier parameter update.
            "mu_max_fact": [float, 1000.0],
            "mu_max": [float, 100000.0],
            "mu_min": [float, 1e-11],
            # Controls how the globalization strategy is applied
            # when the adaptive mu strategy is employed. This
            # controls what quantity is used to control the switch back
            # to monotone mode. Other options are: 'kkt-error', and
            # 'never-monotone-mode' which disables globalization
            "adaptive_mu_globalization": [str, "obj-constr-filter"],
            "adaptive_mu_kkterror_red_iters": [float, 4],
            "adaptive_mu_kkterror_red_fact": [float, 0.9999],
            "filter_margin_fact": [float, 1e-05],
            "filter_max_margin": [float, 1.0],
            "adaptive_mu_restore_previous_iterate": [str, "no"],
            "adaptive_mu_monotone_init_factor": [float, 0.8],
            "adaptive_mu_kkt_norm_type": [str, "2-norm-squared"],
            # Use the mu update strategy: Defaults to
            # Fiacco-McCormick monotone, the other option is
            # 'adaptive'
            "mu_strategy": [str, "monotone"],
            # Select the method used to pick the next mu in
            # the adaptive strategy: Other options: 'loqo':
            # the LOQO adaptive barrier strategy and
            # 'probing': Mehrotra's probing method
            "mu_oracle": [str, "quality-function"],
            "fixed_mu_oracle": [str, "average_compl"],
            # Options for the barrier strategy in IPOPT -
            # these can make a big difference in the
            # performance of the IP algorithm.
            "mu_init": [float, 0.1],
            # Parameter that controls how tightly each barrier
            # problem is solved before the next mu update. A scaled
            # version of: ||KKT||_{infty} < mu*barrier_tol_factor
            "barrier_tol_factor": [float, 10.0],
            # For the monotone strategy, decrease the value of
            # mu by this fixed fraction after each barrrier
            # problem is solved
            "mu_linear_decrease_factor": [float, 0.2],
            # Use the min of the linear decrease factor and
            # mu**(mu_superlinear_decrease_power) for the next
            # barrier value: enables superlinear rates of
            # convergence
            "mu_superlinear_decrease_power": [float, 1.5],
            "mu_allow_fast_monotone_decrease": [str, "yes"],
            "tau_min": [float, 0.99],
            "sigma_max": [float, 100.0],
            "sigma_min": [float, 1e-06],
            "quality_function_norm_type": [str, "2-norm-squared"],
            "quality_function_centrality": [str, "none"],
            "quality_function_balancing_term": [str, "none"],
            "quality_function_max_section_steps": [int, 8],
            "quality_function_section_sigma_tol": [float, 0.01],
            "quality_function_section_qf_tol": [float, 0.0],
            # Line Search.
            "line_search_method": [str, "filter"],
            "alpha_red_factor": [float, 0.5],
            "accept_every_trial_step": [str, "no"],
            "accept_after_max_steps": [int, -1],
            "alpha_for_y": [str, "primal"],
            "alpha_for_y_tol": [float, 10.0],
            "tiny_step_tol": [float, 2.22045e-15],
            "tiny_step_y_tol": [float, 0.01],
            "watchdog_shortened_iter_trigger": [int, 10],
            "watchdog_trial_iter_max": [int, 3],
            "theta_max_fact": [float, 10000.0],
            "theta_min_fact": [float, 0.0001],
            "eta_phi": [float, 1e-08],
            "delta": [float, 1.0],
            "s_phi": [float, 2.3],
            "s_theta": [float, 1.1],
            "gamma_phi": [float, 1e-08],
            "gamma_theta": [float, 1e-05],
            "alpha_min_frac": [float, 0.05],
            # Maximum number of second order correction trial
            # steps at each iteration
            "max_soc": [int, 4],
            "kappa_soc": [float, 0.99],
            "obj_max_inc": [float, 5.0],
            "max_filter_resets": [int, 5],
            "filter_reset_trigger": [int, 5],
            "corrector_type": [str, "none"],
            "skip_corr_if_neg_curv": [str, "yes"],
            "skip_corr_in_monotone_mode": [str, "yes"],
            "corrector_compl_avrg_red_fact": [float, 1.0],
            "nu_init": [float, 1e-06],
            "nu_inc": [float, 0.0001],
            "rho": [float, 0.1],
            "kappa_sigma": [float, 1e10],
            "recalc_y": [str, "no"],
            "recalc_y_feas_tol": [float, 1e-06],
            "slack_move": [float, 1.81899e-12],
            "constraint_violation_norm_type": [str, "1-norm"],
            # Warm Start.
            "warm_start_init_point": [str, "no"],
            "warm_start_same_structure": [str, "no"],
            "warm_start_bound_push": [float, 0.001],
            "warm_start_bound_frac": [float, 0.001],
            "warm_start_slack_bound_push": [float, 0.001],
            "warm_start_slack_bound_frac": [float, 0.001],
            "warm_start_mult_bound_push": [float, 0.001],
            "warm_start_mult_init_max": [float, 1e06],
            "warm_start_entire_iterate": [float, "no"],
            # Linear Solver.
            "linear_solver": [str, "mumps"],
            "linear_system_scaling": [str, "none"],  # Had been "mc19", but not always available.
            "linear_scaling_on_demand": [str, "yes"],
            # Step Calculation.
            # Use Mehrotra's predictor-corrector algorithm -
            # warning: no globalization
            "mehrotra_algorithm": [str, "no"],
            "fast_step_computation": [str, "no"],
            "min_refinement_steps": [int, 1],
            "max_refinement_steps": [int, 10],
            "residual_ratio_max": [float, 1e-10],
            "residual_ratio_singular": [float, 1e-05],
            "residual_improvement_factor": [float, 1.0],
            "neg_curv_test_tol": [float, 1.0],
            "max_hessian_perturbation": [float, 1e20],
            "min_hessian_perturbation": [float, 1e-20],
            "perturb_inc_fact_first": [float, 100.0],
            "perturb_inc_fact": [float, 8.0],
            "perturb_dec_fact": [float, 0.333333],
            "first_hessian_perturbation": [float, 0.0001],
            "jacobian_regularization_value": [float, 1e-08],
            "jacobian_regularization_exponent": [float, 0.25],
            "perturb_always_cd": [str, "no"],
            # Restoration Phase.
            "expect_infeasible_problem": [str, "no"],
            "expect_infeasible_problem_ctol": [float, 0.001],
            "expect_infeasible_problem_ytol": [float, 1e08],
            "start_with_resto": [str, "no"],
            "soft_resto_pderror_reduction_factor": [float, 0.9999],
            "max_soft_resto_iters": [int, 10],
            "required_infeasibility_reduction": [float, 0.9],
            "max_resto_iter": [int, 3000000],
            "evaluate_orig_obj_at_resto_trial": [str, "yes"],
            "resto_penalty_parameter": [float, 1000.0],
            "resto_proximity_weight": [float, 1.0],
            "bound_mult_reset_threshold": [float, 1000.0],
            "constr_mult_reset_threshold": [float, 0.0],
            "resto_failure_feasibility_threshold": [float, 0.0],
            # Derivative Checker.
            # Derivative Testing options
            # none, first-order, second-order, only-second-order
            "derivative_test": [str, "none"],
            "derivative_test_first_index": [int, -2],
            "derivative_test_perturbation": [float, 1e-08],
            "derivative_test_tol": [float, 1e-4],
            "derivative_test_print_all": [str, "no"],
            "jacobian_approximation": [str, "exact"],
            "findiff_perturbation": [float, 1e-07],
            "point_perturbation_radius": [float, 10.0],
            # Hessian Approximation.
            "limited_memory_aug_solver": [str, "sherman-morrison"],
            "limited_memory_max_history": [int, 6],
            "limited_memory_update_type": [str, "bfgs"],
            "limited_memory_initialization": [str, "scalar1"],
            "limited_memory_init_val": [float, 1.0],
            "limited_memory_init_val_max": [float, 1e08],
            "limited_memory_init_val_min": [float, 1e-08],
            "limited_memory_max_skipping": [int, 2],
            "limited_memory_special_for_resto": [str, "no"],
            "hessian_approximation": [str, "limited-memory"],
            "hessian_approximation_space": [str, "nonlinear-variables"],
            # MA27 Linear Solver.
            "ma27_pivtol": [float, 1e-08],
            "ma27_pivtolmax": [float, 0.0001],
            "ma27_liw_init_factor": [float, 5.0],
            "ma27_la_init_factor": [float, 5.0],
            "ma27_meminc_factor": [float, 10.0],
            "ma27_skip_inertia_check": [str, "no"],
            "ma27_ignore_singularity": [str, "no"],
            # # MA57 Linear Solver.
            # 'ma57_pivtol' : [float, 1e-08],
            # 'ma57_pivtolmax' : [float, 0.0001],
            # 'ma57_pre_alloc' : [float, 1.05],
            # 'ma57_pivot_order' : [int, 5],
            # 'ma57_automatic_scaling' : [str, "yes"],  # ipopt default is "no".
            # 'ma57_block_size' : [int, 16],
            # 'ma57_node_amalgamation' : [int, 16],
            # 'ma57_small_pivot_flag' : [int, 0],
            # # Paridiso Linear Solver.
            # 'pardiso_matching_strategy' : [str, "complete+2x2"],
            # 'pardiso_redo_symbolic_fact_only_if_inertia_wrong' : [str, "no"],
            # 'pardiso_repeated_perturbation_means_singular' : [str, "no"],
            # 'pardiso_msglvl' : [int, 0],
            # 'pardiso_skip_inertia_check' : [str, "no"],
            # 'pardiso_max_iter' : [int, 500],
            # 'pardiso_iter_relative_tol' : [float, 1e-06],
            # 'pardiso_iter_coarse_size' : [int, 5000],
            # 'pardiso_iter_max_levels' : [int, 10000],
            # 'pardiso_iter_dropping_factor' : [float, 0.5],
            # 'pardiso_iter_dropping_schur' : [float, 0.1],
            # 'pardiso_iter_max_row_fill' : [int, 10000000],
            # 'pardiso_iter_inverse_norm_factor' : [float, 5e+06],
            # 'pardiso_iterative' : [str, "no"],
            # 'pardiso_max_droptol_corrections' : [int, 4],
            # # Mumps Linear Solver.
            # 'mumps_pivtol' : [float, 1e-06],
            # 'mumps_pivtolmax' : [float, 0.1],
            # 'mumps_mem_percent' : [int, 1000],
            # 'mumps_permuting_scaling' : [int, 7],
            # 'mumps_pivot_order' : [int, 7],
            # 'mumps_scaling' : [int, 77],
            # 'mumps_dep_tol' : [float, -1.0],
            # # MA28 Linear Solver.
            # 'ma28_pivtol' : [float, 0.01],
            # # Uncategorized.
            # 'warm_start_target_mu' : [float, 0.0]
        }

        self.informs = {
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

        if pyipoptcore is None:
            if raiseError:
                raise Error("There was an error importing the compiled IPOPT module")

        self.set_options = []
        Optimizer.__init__(self, name, category, self.defOpts, self.informs, *args, **kwargs)

        # IPOPT needs Jacobians in coo format
        self.jacType = "coo"

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

        self.callCounter = 0
        self.storeSens = storeSens

        if len(optProb.constraints) == 0:
            # If the user *actually* has an unconstrained problem,
            # snopt sort of chokes with that....it has to have at
            # least one constraint. So we will add one
            # automatically here:
            self.unconstrained = True
            optProb.dummyConstraint = True

        # Save the optimization problem and finalize constraint
        # Jacobian, in general can only do on root proc
        self.optProb = optProb
        self.optProb.finalizeDesignVariables()
        self.optProb.finalizeConstraints()
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
            blc = np.array([-1e20])
            buc = np.array([1e20])
            ncon = 1

        jac = convertToCOO(jac)  # Conver to coo format for IPOPT

        # We make a split here: If the rank is zero we setup the
        # problem and run IPOPT, otherwise we go to the waiting loop:

        if self.optProb.comm.rank == 0:

            # Now what we need for IPOPT is precisely the .row and
            # .col attributes of the fullJacobian array
            matStruct = (jac["coo"][IROW].copy().astype("int_"), jac["coo"][ICOL].copy().astype("int_"))

            # Set history/hotstart
            self._setHistory(storeHistory, hotStart)

            # Define the 4 call back functions that ipopt needs:
            def eval_f(x, user_data=None):
                fobj, fail = self._masterFunc(x, ["fobj"])
                return fobj

            def eval_g(x, user_data=None):
                fcon, fail = self._masterFunc(x, ["fcon"])
                return fcon.copy()

            def eval_grad_f(x, user_data=None):
                gobj, fail = self._masterFunc(x, ["gobj"])
                return gobj.copy()

            def eval_jac_g(x, flag, user_data=None):
                if flag:
                    return copy.deepcopy(matStruct)
                else:
                    gcon, fail = self._masterFunc(x, ["gcon"])
                    return gcon.copy()

            timeA = time.time()
            nnzj = len(matStruct[0])
            nnzh = 0

            nlp = pyipoptcore.create(
                len(xs), blx, bux, ncon, blc, buc, nnzj, nnzh, eval_f, eval_grad_f, eval_g, eval_jac_g
            )

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
        set all of the the options in self.set_options in the ipopt instance nlp
        """
        # Set Options from the local options dictionary
        # ---------------------------------------------

        for key in self.options:
            if key != "defaults":
                name = key
                value = self.getOption(key)

                if isinstance(value, str):
                    nlp.str_option(name, value)
                elif isinstance(value, float):
                    nlp.num_option(name, value)
                elif isinstance(value, int):
                    nlp.int_option(name, value)
                else:
                    print("invalid option type", type(value))

    def _on_setOption(self, name, value):
        """
        Set Optimizer Option Value (Optimizer Specific Routine)

        Documentation last updated:  May. 07, 2008 - Ruben E. Perez
        """

        self.set_options.append([name, value])

    def _on_getOption(self, name):
        """
        Routine to be implemented by optimizer
        """

        pass

    def _on_getInform(self, info):
        """
        Routine to be implemented by optimizer
        """

        pass


# ==============================================================================
# IPOPT Optimizer Test
# ==============================================================================
if __name__ == "__main__":

    ipopt = IPOPT()
    print(ipopt)
