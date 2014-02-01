#from __future__ import absolute_import
#/bin/env python
'''
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
'''
# =============================================================================
# IPOPT Library
# =============================================================================
#from . import pyipoptcore
import pyipoptcore
# try:
#     import pyipoptcore
# except:
#     raise ImportError('IPOPT shared library failed to import')

# =============================================================================
# Standard Python modules
# =============================================================================
import os
import sys
import copy
import time
import types
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
from ..pyOpt_gradient import Gradient
from ..pyOpt_solution import Solution
from ..pyOpt_error import Error
# =============================================================================
# Misc Definitions
# =============================================================================
inf = 1e20  # define a value for infinity

# Try to import mpi4py and determine rank
try: 
    from mpi4py import MPI
    rank = MPI.COMM_WORLD.rank
except:
    rank = 0
    MPI = None
# end try

# =============================================================================
# IPOPT Optimizer Class
# =============================================================================
class IPOPT(Optimizer):
    '''
    IPOPT Optimizer Class - Inherited from Optimizer Abstract Class
    '''
    
    def __init__(self, *args, **kwargs):
        '''
        IPOPT Optimizer Class Initialization
        '''
        
        name = 'IPOPT'
        category = 'Local Optimizer'
        def_opts = {# Don't have any yet....
            }
        informs = { # Don't have any of these yet either..
            }

        self.set_options = []
        Optimizer.__init__(self, name, category, def_opts, informs, *args, **kwargs)

        # IPOPT needs jacobians in coo format
        self.jacType = 'coo'

        # Constrained until we know otherwise :-)
        self.unconstrained = False

    def __call__(self, optProb, sens=None, sensStep=None, sensMode=None,
                  storeHistory=None, hotStart=None, 
                  coldStart=None, timeLimit=None, comm=None):
        '''
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

        timeLimit : number
            Number of seconds to run the optimization before a
            terminate flag is given to the optimizer and a "clean"
            exit is performed.

        comm : MPI Intra communicator
            Specifiy a MPI comm to use. Default is None. If mpi4py is not
            available, the serial mode will still work. if mpi4py *is*
            available, comm defaluts to MPI.COMM_WORLD. 
            '''
        
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
        if self.optProb.nlCon > 0:
            self.appendLinearConstraints

        # Setup initial cache values
        self._setInitialCacheValues()

        # Next we determine what to what to do about
        # derivatives. We must hvae a function or we use FD or CS:
        if isinstance(sens, types.FunctionType):
            # We have function handle for gradients! Excellent!
            self.sens = sens
        elif sens.lower() in ['fd','cs']:
            # Create the gradient class that will operate just like if
            # the user supplied fucntion
            self.sens = Gradient(optProb, sens.lower(), sensStep,
                                 sensMode, comm)
        else:
            raise Error('Unknown value given for sens. Must be None, \'FD\', \
            \'CS\' or a python function handle')
        # end if
                
        # We make a split here: If the rank is zero we setup the
        # problem and run IPOPT, otherwise we go to the waiting loop:

        if rank == 0:

            # Get the variable names and variable bounds
            # ------------------------------------------
            blx = []
            bux = []
            xs = []
            for dvSet in self.optProb.variables.keys():
                for dvGroup in self.optProb.variables[dvSet]:
                    for var in self.optProb.variables[dvSet][dvGroup]:
                        if var.type == 'c':
                            blx.append(var.lower)
                            bux.append(var.upper)
                            xs.append(var.value)

            blx = numpy.array(blx)
            bux = numpy.array(bux)
            xs = numpy.array(xs)

            # Constraints Handling -- Put the nonlinear constraints
            # first which will make it easier to tack on the linear
            # ones at the end.
            blc = []
            buc = []
                
            for key in self.optProb.constraints.keys():
                if not self.optProb.constraints[key].linear:
                    blc.extend(self.optProb.constraints[key].lower)
                    buc.extend(self.optProb.constraints[key].upper)
         
            for key in self.optProb.constraints.keys():
                if self.optProb.constraints[key].linear:
                    blc.extend(self.optProb.constraints[key].lower)
                    buc.extend(self.optProb.constraints[key].upper)
            if self.unconstrained:
                blc.append(-inf)
                buc.append(inf)

            ncon = len(blc)
            blc = numpy.array(blc)
            buc = numpy.array(buc)

            # Before we start, we assemble the full jacobian, convert
            # to COO format, and store the format since we will need
            # that on the first constraint jacobian call back. 
            # Determine the sparsity structure of the full jacobian
            # -----------------------------------------------------
            # Get nonlinear part:
            gcon = {}
            for iCon in self.optProb.constraints:
                con = self.optProb.constraints[iCon]
                if not con.linear:
                    gcon[iCon] = con.jac

            fullJacobian = self.optProb.processConstraintJacobian(
                gcon, linearFlag=False)

            # If we have linear constraints those are already done actually.
            if self.optProb.linearJacobian is not None:
                fullJacobian = sparse.vstack([fullJacobian,
                                              self.optProb.linearJacobian])

            # Now what we need for IPOPT is precisely the .row and
            # .col attributes of the fullJacobian array
            matStruct = (fullJacobian.row.copy().astype('int64'), 
                         fullJacobian.col.copy().astype('int64'))

            self._setHistory(storeHistory)
            self._hotStart(storeHistory, hotStart)

            # Define the 4 call back functions that ipopt needs:
            def eval_f(x, user_data=None):
                fobj, fail = self.masterFunc(x, ['fobj'])
                return fobj

            def eval_g(x, user_data = None):
                fcon, fail = self.masterFunc(x, ['fcon'])
                return fcon.copy()

            def eval_grad_f(x, user_data= None):
                gobj, fail = self.masterFunc(x, ['gobj'])
                return gobj.copy()

            def eval_jac_g(x, flag, user_data = None):
                if flag:
                    print matStruct
                    return copy.deepcopy(matStruct)
                else:
                    gcon, fail = self.masterFunc(x, ['gcon'])
                    return gcon.copy()
            
            timeA = time.time()
            nnzj = len(matStruct[0])
            nnzh = 0
            nlp = pyipoptcore.create(len(xs), blx, bux, ncon, blc, buc, nnzj, nnzh, 
                                          eval_f, eval_grad_f, eval_g, eval_jac_g) 
            nlp.num_option('tol',1e-5)
            nlp.str_option('hessian_approximation','limited-memory')
            nlp.int_option('limited_memory_max_history',10)
            nlp.str_option('derivative_test','first-order')
            nlp.int_option('max_iter',10)
            x, zl, zu, constraint_multipliers, obj, status = nlp.solve(xs)
            nlp.close()
            optTime = time.time()-timeA
            
            # Store Results
            sol_inform = {}
            # sol_inform['value'] = inform
            # sol_inform['text'] = self.informs[inform[0]]

            # Create the optimization solution
            sol = Solution(self.optProb, optTime, 1, sol_inform)

            # Now set the x-values:
            i = 0
            for dvSet in sol.variables.keys():
                for dvGroup in sol.variables[dvSet]:
                    for var in sol.variables[dvSet][dvGroup]:
                        var.value = xs[i]
                        i += 1
            sol.fStar = obj

            if MPI:
                # Broadcast a -1 to indcate IPOPT has finished
                MPI.COMM_WORLD.bcast(-1, root=0)

        else:
            self.waitLoop()
            sol = None
        # end if

        # Communicate the solution -- We are back to the point where
        # all processors are back together, so a standard bcast is
        # fine.
        if MPI:
            sol = MPI.COMM_WORLD.bcast(sol)

        return  sol


    def _on_setOption(self, name, value):
        '''
        Set Optimizer Option Value (Optimizer Specific Routine)
        
        Documentation last updated:  May. 07, 2008 - Ruben E. Perez
        '''
        
        self.set_options.append([name,value])
        
    def _on_getOption(self, name):
        '''
        Get Optimizer Option Value (Optimizer Specific Routine)
        
        Documentation last updated:  May. 07, 2008 - Ruben E. Perez
        '''
        
        pass
        
    def _on_getInform(self, infocode):
        '''
        Get Optimizer Result Information (Optimizer Specific Routine)
        
        Keyword arguments:
        -----------------
        id -> STRING: Option Name
        
        Documentation last updated:  May. 07, 2008 - Ruben E. Perez
        '''
        
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
        '''
        Flush the Output Files (Optimizer Specific Routine)
        
        Documentation last updated:  August. 09, 2009 - Ruben E. Perez
        '''
        
        pass
 
#==============================================================================
# IPOPT Optimizer Test
#==============================================================================
if __name__ == '__main__':
    
    ipopt = IPOPT()
    print ipopt
