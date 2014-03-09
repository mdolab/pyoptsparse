#/bin/env python
"""
pyIPOPT - A python wrapper to the core IPOPT compiled module.

Copyright (c) 2013-2014 by Dr. Gaetan Kenway
All rights reserved.

Tested on:
---------
Linux with intel

Developers:
-----------
- Dr. Gaetan Kenway (GKK)
- Dr. Graeme Kennedy (GJK)
History
-------
    v. 0.1    - Initial Wrapper Creation
"""
# =============================================================================
# IPOPT Library
# =============================================================================

try:
    from . import pyipoptcore
except:
    pyipoptcore = None

# =============================================================================
# standard Python modules
# =============================================================================
import os
import copy
import time
# =============================================================================
# External Python modules
# =============================================================================
import numpy
import shelve
from scipy import sparse
# # =============================================================================
# # Extension modules
# # =============================================================================
from ..pyOpt_optimizer import Optimizer
from ..pyOpt_history import History
from ..pyOpt_error import Error
# =============================================================================
# IPOPT Optimizer Class
# =============================================================================
class IPOPT(Optimizer):
    """
    IPOPT Optimizer Class - Inherited from Optimizer Abstract Class
    """

    def __init__(self, *args, **kwargs):
        """
        IPOPT Optimizer Class Initialization
        """

        name = 'IPOPT'
        category = 'Local Optimizer'
        def_opts = {'tol':[float,1e-6],
                    'hessian_approximation':[str,'limited-memory'],
                    'limited_memory_max_history':[int,10],
                    'max_iter':[int,100],
                    # print options
                    'print_level':[int, 0], # Output verbosity level. '0-12'
                    'print_user_options':[str,'no'], #yes or no, Print all options set by the user.
                    'print_options_documentation':[str,'no'],#yes or no,Switch to print all algorithmic options.
                    'print_frequency_iter':[int,1],# Determines at which iteration frequency the summarizing iteration output line should be printed.
                    'print_frequency_time':[float, 0.0],# Determines at which time frequency the summarizing iteration output line should be printed. 
                    'output_file':[str,'IPOPT.out'],
                    'file_print_level':[int,5],#Verbosity level for output file. '0-12'
                    'option_file_name':[str,'IPOPT_options.opt'],
                    'print_info_string':[str,'no'],#yes or no.Enables printing of additional info string at end of iteration output.
                    'inf_pr_output':[str,'original'],#Determines what value is printed in the "inf_pr" output column. 'internal' or 'original'
                    'print_timing_statistics':[str,'no'],#yes or no
                    # Derivative Testing options
                    'derivative_test':[str,'none'], # none,first-order,second-order,only-second-order
                    'derivative_test_perturbation':[float,1e-8],
                    'derivative_test_tol':[float,1e-4],
                    'derivative_test_print_all':[str,'no'],#yes,no
                    'derivative_test_first_index':[int,-2],
                    'point_perturbation_radius':[float, 10.0],
                    # Line search
                    'max_soc':[int, 4], #Maximum numbero fsecond order correction trial steps at each iteration
                    'watchdog_shortened_iter_trigger':[int, 10],
                    'watchdog_trial_iter_max':[int,3],
                    'accept_every_trial_step':[str, 'no'],
                    'corrector_type':[str,'none'],
                    'mu_init':[float, 0.1],
                    'mu_strategy':[str,'monotone'],
                    'start_with_resto':[str,'no'],
                    'required_infeasibility_reduction':[float,0.9],
                    'expect_infeasible_problem':[str,'no'],
                    }
        informs = { # Don't have any of these yet either..
            }

        if pyipoptcore is None:
            raise Error('There was an error importing the compiled \
                        IPOPT module')

        self.set_options = []
        Optimizer.__init__(self, name, category, def_opts, informs, *args, **kwargs)

        # IPOPT needs jacobians in coo format
        self.jacType = 'coo'

    def __call__(self, optProb, sens=None, sensStep=None, sensMode=None,
                  storeHistory=None, hotStart=None,
                  coldStart=None):
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

        coldStart : str
            Filename of the history file to use for "cold"
            restart. Here, the only requirment is that the number of
            design variables (and their order) are the same. Use this
            method if any of the optimization parameters have changed.
            """
       
        self.callCounter = 0

        if len(optProb.constraints) == 0:
            # If the user *actually* has an unconstrained problem,
            # snopt sort of chokes with that....it has to have at
            # least one constraint. So we will add one
            # automatically here:
            self.unconstrained = True
            optProb.dummyConstraint = True

        # Save the optimization problem and finialize constraint
        # jacobian, in general can only do on root proc
        self.optProb = optProb
        self.optProb.finalizeDesignVariables()
        self.optProb.finalizeConstraints()
        self._setInitialCacheValues()
        blx, bux, xs = self._assembleContinuousVariables()
        self._setSens(sens, sensStep, sensMode)
        ff = self._assembleObjective()

        # Determine the sparsity structure of the full jacobian
        # -----------------------------------------------------

        # Gather dummy data and process jacobian:
        gcon = {}
        for iCon in self.optProb.constraints:
            gcon[iCon] = self.optProb.constraints[iCon].jac

        jac = self.optProb.processConstraintJacobian(gcon)

        if self.optProb.nCon > 0:
            # We need to reorder this full jacobian...so get ordering:
            indices, blc, buc, fact = self.optProb.getOrdering(
                ['ne','ni','le','li'], oneSided=False)
            self.optProb.jacIndices = indices
            self.optProb.fact = fact
            self.optProb.offset = numpy.zeros(len(indices))
            ncon = len(indices)
            jac = jac[indices, :] # Does reordering
            jac = fact*jac # Perform logical scaling
        else:
            blc = numpy.array([-1e20])
            buc = numpy.array([1e20])
            ncon = 1
        jac = jac.tocoo() # Conver to coo format for IPOPT

        # We make a split here: If the rank is zero we setup the
        # problem and run IPOPT, otherwise we go to the waiting loop:

        if self.optProb.comm.rank == 0:

            # Now what we need for IPOPT is precisely the .row and
            # .col attributes of the fullJacobian array
            matStruct = (jac.row.copy().astype('int_'),
                         jac.col.copy().astype('int_'))
            self._setHistory(storeHistory)
            self._hotStart(storeHistory, hotStart)

            # Define the 4 call back functions that ipopt needs:
            def eval_f(x, user_data=None):
                fobj, fail = self._masterFunc(x, ['fobj'])
                return fobj

            def eval_g(x, user_data = None):
                fcon, fail = self._masterFunc(x, ['fcon'])
                return fcon.copy()

            def eval_grad_f(x, user_data= None):
                gobj, fail = self._masterFunc(x, ['gobj'])
                return gobj.copy()

            def eval_jac_g(x, flag, user_data = None):
                if flag:
                    return copy.deepcopy(matStruct)
                else:
                    gcon, fail = self._masterFunc(x, ['gcon'])
                    return gcon.copy()

            timeA = time.time()
            nnzj = len(matStruct[0])
            nnzh = 0

            nlp = pyipoptcore.create(len(xs), blx, bux, ncon, blc, buc, nnzj, nnzh,
                                     eval_f, eval_grad_f, eval_g, eval_jac_g)

            self._set_ipopt_options(nlp)
            x, zl, zu, constraint_multipliers, obj, status = nlp.solve(xs)
            nlp.close()
            optTime = time.time()-timeA

            if self.storeHistory:
                self.hist.close()

            # Store Results
            sol_inform = {}
            # sol_inform['value'] = inform
            # sol_inform['text'] = self.informs[inform[0]]

            # Create the optimization solution
            sol = self._createSolution(optTime, sol_inform, obj, x)

            # Indicate solution finished
            self.optProb.comm.bcast(-1, root=0)
        else:
            self._waitLoop()
            sol = None

        # Communication solution and return
        sol = self._communicateSolution(sol)

        return  sol

    def _set_ipopt_options(self, nlp):
        """
        set all of the the options in self.set_options in the ipopt instance nlp
        """
        # Set Options from the local options dictionary
        # ---------------------------------------------

        for key in self.options:
            if key != 'defaults':
                name = key
                value = self.getOption(key)

                if isinstance(value, str):
                    nlp.str_option(name,value)
                elif isinstance(value, float):
                    nlp.num_option(name,value)
                elif isinstance(value, int):
                    nlp.int_option(name,value)
                else:
                    print 'invalid option type',type(value)


    def _on_setOption(self, name, value):
        """
        Set Optimizer Option Value (Optimizer Specific Routine)

        Documentation last updated:  May. 07, 2008 - Ruben E. Perez
        """

        self.set_options.append([name,value])

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
        mjr_code = (infocode[0]/10)*10
        mnr_code = infocode[0] - 10*mjr_code
        try:
            inform_text = self.informs[mjr_code]
        except:
            inform_text = 'Unknown Exit Status'
        # end try

        return inform_text

    def _on_flushFiles(self):
        """
        Flush the Output Files (Optimizer Specific Routine)

        Documentation last updated:  August. 09, 2009 - Ruben E. Perez
        """

        pass

#==============================================================================
# IPOPT Optimizer Test
#==============================================================================
if __name__ == '__main__':

    ipopt = IPOPT()
    print ipopt
