#/bin/env python
"""
pyNLPY_AUGLAG - A wrapper for the matrix-free augmented Lagrangian algorithm 
from the NLPy package. This optimizer can exploit sparsity in the constraint 
blocks.

Some assumptions to get this wrapper working:
- the NLPy package has been downloaded from Github and installed locally
- only one version of the optimizer is available: LSR1 Hessian and adjoint 
Broyden A Jacobian approximation (though this is the best in our experience)

Copyright (c) 2014 by Andrew Lambe
All rights reserved.

Tested on:
---------
Linux with intel

Developers:
-----------
- Andrew Lambe (ABL)
History
-------
    v. 0.1    - Initial Wrapper Creation
"""
from __future__ import absolute_import
from __future__ import print_function
# =============================================================================
# NLPy Library
# =============================================================================
from nlpy.model.mfnlp import MFModel
from nlpy.optimize.solvers.sbmin import SBMINTotalLqnFramework
from nlpy.optimize.solvers.auglag2 import AugmentedLagrangianTotalLsr1AdjBroyAFramework
# =============================================================================
# Standard Python modules
# =============================================================================
import os
import time
import copy
# =============================================================================
# External Python modules
# =============================================================================
import numpy
eps = numpy.finfo(1.0).eps
import logging
from mpi4py import MPI
# ===========================================================================
# Extension modules
# ===========================================================================
from ..pyOpt_optimizer import Optimizer
from ..pyOpt_gradient import Gradient
from ..pyOpt_history import History
# from ..pyOpt_error import Error
# =============================================================================
# NLPy Optimizer Class
# =============================================================================
class NLPY_AUGLAG(Optimizer):
    """
    NLPY_AUGLAG Optimizer Class - Inherited from Optimizer Abstract Class
    """

    def __init__(self, *args, **kwargs):
        """
        NLPY_AUGLAG Optimizer Class Initialization
        """

        name = 'NLPY_AUGLAG'
        category = 'Local Optimizer'
        # Default options go here - many to be added
        defOpts = {
        'Logger Name':[str,'nlpy_logging'],
        'Prefix':[str,'./'],
        'Save Current Point':[bool,False],
        'Absolute Optimality Tolerance':[float,1.0e-6],
        'Relative Optimality Tolerance':[float,1.0e-6],
        'Absolute Feasibility Tolerance':[float,1.0e-6],
        'Relative Feasibility Tolerance':[float,1.0e-6],
        'Number of Quasi-Newton Pairs':[int,5],
        'Use N-Y Backtracking':[bool,True],
        'Use Magical Steps':[bool,True],
        }
        # Inform/Status codes go here
        informs = {
        0: 'Successfully converged',
        -1: 'Maximum number of iterations reached',
        -2: 'Problem appears to be infeasible',
        -3: 'Solver stopped on user request',
        }
        self.set_options = []
        Optimizer.__init__(self, name, category, defOpts, informs, *args, **kwargs)


    def __call__(self, optProb, sens=None, sensStep=None, sensMode=None,
                  storeHistory=None, hotStart=None,
                  coldStart=None, timeLimit=None):
        """
        This is the main routine used to call the optimizer.

        Parameters
        ----------
        optProb : Optimization or Solution class instance
            This is the complete description of the optimization problem
            to be solved by the optimizer

        sens : str or python Function or list of functions.
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
        """

        self.callCounter = 0

        # NLPy *can* handle unconstrained problems cleanly
        # Check if we still need this later
        # if len(optProb.constraints) == 0:
        #     self.unconstrained = True
        #     optProb.dummyConstraint = True

        # Save the optimization problem and finialize constraint
        # jacobian, in general can only do on root proc
        self.optProb = optProb
        self.optProb.finalizeDesignVariables()
        self.optProb.finalizeConstraints()

        self._setInitialCacheValues()
        blx, bux, xs = self._assembleContinuousVariables()

        # Need to fix the "set sensitivity" function
        self._setSens(sens, sensStep, sensMode)
        ff = self._assembleObjective()

        oneSided = False

        # Need both constraint orderings and variable orderings here
        # to correctly compute matrix-vector products

        # gcon = {}
        # for iCon in self.optProb.constraints:
        #     gcon[iCon] = self.optProb.constraints[iCon].jac

        # jac = self.optProb.processConstraintJacobian(gcon)

        if self.optProb.nCon > 0:
            # We need to reorder this full jacobian...so get ordering:
            indices, blc, buc, fact = self.optProb.getOrdering(
                ['ne','ni','le','li'], oneSided=False)
            self.optProb.jacIndices = indices
            self.optProb.fact = fact
            self.optProb.offset = numpy.zeros(len(indices))
            ncon = len(indices)
            # jac = jac[indices, :] # Does reordering
            # jac = fact*jac # Perform logical scaling
        else:
            blc = numpy.array([])
            buc = numpy.array([])
            # ncon = 1

        # Since this algorithm exploits parallel computing, define the 
        # callbacks on all processors
        def obj(x):
            fobj, fail = self._masterFunc(x, ['fobj'])
            return fobj

        def cons(x):
            fcon, fail = self._masterFunc(x, ['fcon'])
            return fcon.copy()

        # Gradient callbacks are specialized for the matrix-free case
        if self.matrix_free == True:

            # Note the different tags for calling matrix-free and 
            # standard matrix-vector product functions

            # Also note the different 'masterFunc' call for the objective 
            # gradient. This is necessary as the structure of the callback 
            # functions is fundamentally different.

            def grad(x):
                gobj, fail = self._masterFunc3(x, None, ['gobj'])
                return gobj.copy()

            def jprod(x,p,sparse_only=False):
                q, fail = self._masterFunc3(x, p, ['jac_prod'], sparse_only=sparse_only)
                return q.copy()

            def jtprod(x,q,sparse_only=False):
                p, fail = self._masterFunc3(x, q, ['jac_T_prod'], sparse_only=sparse_only)
                return p.copy()

        else:

            def grad(x):
                gobj, fail = self._masterFunc(x, ['gobj'])
                return gobj.copy()

            def jprod(x,p,sparse_only=False):
                q, fail = self._masterFunc3(x, p, ['gcon_prod'], sparse_only=sparse_only)
                return q.copy()

            def jtprod(x,q,sparse_only=False):
                p, fail = self._masterFunc3(x, q, ['gcon_T_prod'], sparse_only=sparse_only)
                return p.copy()

        # end if

        # Step 2: Set up the optimizer with all the necessary functions
        nlpy_problem = MFModel(n=self.optProb.ndvs,m=self.optProb.nCon,name=optProb.name,x0=xs,
            Lvar=blx, Uvar=bux, Lcon=blc, Ucon=buc)
        nlpy_problem.obj = obj
        nlpy_problem.cons = cons
        nlpy_problem.grad = grad
        nlpy_problem.jprod = jprod
        nlpy_problem.jtprod = jtprod

        # Set up the loggers on one proc only
        if self.optProb.comm.rank == 0:
            lprefix = self.options['Prefix'][1]
            lname = self.options['Logger Name'][1]

            fmt = logging.Formatter('%(name)-15s %(levelname)-8s %(message)s')
            hndlr = logging.FileHandler(lprefix+lname+'_more.log',mode='w')
            hndlr.setLevel(logging.DEBUG)
            hndlr.setFormatter(fmt)
            hndlr2 = logging.FileHandler(lprefix+lname+'.log',mode='w')
            hndlr2.setLevel(logging.INFO)
            hndlr2.setFormatter(fmt)

            # Configure auglag logger.
            auglaglogger = logging.getLogger('nlpy.auglag')
            auglaglogger.setLevel(logging.DEBUG)
            auglaglogger.addHandler(hndlr)
            auglaglogger.addHandler(hndlr2)
            auglaglogger.propagate = False

            # Configure sbmin logger.
            sbminlogger = logging.getLogger('nlpy.sbmin')
            sbminlogger.setLevel(logging.DEBUG)
            sbminlogger.addHandler(hndlr)
            sbminlogger.addHandler(hndlr2)
            sbminlogger.propagate = False

            # Configure bqp logger.
            bqplogger = logging.getLogger('nlpy.bqp')
            bqplogger.setLevel(logging.INFO)
            bqplogger.addHandler(hndlr)
            bqplogger.propagate = False

            # pyOpt History logging for hot starts
            self._setHistory(storeHistory)
            self._hotStart(storeHistory, hotStart)
        # end if

        # Step 3: Pass options and solve the problem
        timeA = time.time()

        # Also need to pass the number of dense constraints
        # Assume the dense block is listed first in the problem definition
        solver = AugmentedLagrangianTotalLsr1AdjBroyAFramework(nlpy_problem, 
            SBMINTotalLqnFramework, 
            omega_abs=self.options['Absolute Optimality Tolerance'][1], 
            eta_abs=self.options['Absolute Feasibility Tolerance'][1], 
            omega_rel=self.options['Relative Optimality Tolerance'][1],
            eta_rel=self.options['Relative Feasibility Tolerance'][1],
            qn_pairs=self.options['Number of Quasi-Newton Pairs'][1],
            data_prefix=self.options['Prefix'][1],
            save_data=self.options['Save Current Point'][1])
            # sparse_index=struct_opt.num_dense_con,
            # data_prefix=prefix, save_data=False)

        solver.solve(ny=self.options['Use N-Y Backtracking'], magic_steps_agg=self.options['Use Magical Steps'])

        # Step 4: Collect and return solution
        optTime = time.time() - timeA

        if self.storeHistory:
            self.hist.close()

        sol_inform = {}
        sol_inform['value'] = solver.status
        sol_inform['text'] = self.informs[solver.status]
        sol = self._createSolution(optTime, sol_inform, solver.f)

        return sol


    def _setSens(self, sens, sensStep, sensMode):
        """
        For the matrix-free approach, the sens argument is actually a list 
        of three separate functions. The order of these functions must be 
        [obj_grad, jac_prod, jac_t_prod].

        Otherwise, this function is identical to that of the base class.
        """

        self.matrix_free = False
        if sens is None:
            raise Error('\'None\' value given for sens. Must be one \
of \'FD\' or \'CS\' or a user supplied function or group of functions.')
        elif hasattr(sens, 'append'):
            # A list of functions has been provided
            self.sens = sens
            self.matrix_free = True
        elif hasattr(sens, '__call__'):
            # A single function has been provided, old-style sensitivities
            self.sens = sens
        elif sens.lower() in ['fd', 'cs']:
            # Create the gradient class that will operate just like if
            # the user supplied fucntion
            self.sens = Gradient(self.optProb, sens.lower(), sensStep,
                                 sensMode, self.optProb.comm)
        else:
            raise Error('Unknown value given for sens. Must be None, \'FD\', \
            \'CS\', a python function handle, or a list of handles')


    def _masterFunc3(self, x, invec, evaluate, writeHist=True, sparse_only=False):
        """
        A shell function for the matrix-free case, called on all processors.

        ** Right now, we assume that the 'evaluate' list has only one element **
        """

        xscaled = self.optProb.invXScale.dot(x)
        xuser = self.optProb.processX(xscaled)

        masterFail = False

        # History storage is a little different for the matrix-free case.
        # We store the sequence of input and output vectors called for in 
        # the optimizer. Storage of the objective function gradient is 
        # unchanged.

        # Set basic parameters in history
        hist = {'x': x, 'xuser': xuser, 'xscaled': xscaled}
        returns = []

        # Evaluate the gradient of the objective function only
        if 'gobj' in evaluate:
            if numpy.linalg.norm(x-self.cache['x']) > eps:
                # Previous evaluated point is *different* than the
                # point requested for the derivative. Recusively call
                # the routine with ['fobj', and 'fcon']
                self._masterFunc2(x, ['fobj', 'fcon'], writeHist=False)

            if self.cache['gobj'] is None:
                # Only the objective function is cached in this version
                timeA = time.time()
                gobj, fail = self.sens[0](xuser)
                self.userSensTime += time.time()-timeA
                self.userSensCalls += 1
                # User values are stored immediately
                self.cache['gobj_user'] = copy.deepcopy(gobj)

                # Process objective gradient for optimizer
                gobj = self.optProb.processObjectiveGradient(gobj)

                # Set the cache values
                self.cache['gobj'] = gobj.copy()

                # Update fail flag
                masterFail = masterFail or fail

            # gobj is now in the cache
            returns.append(self.cache['gobj'])
            hist['gobj'] = self.cache['gobj_user']

        # Evaluate the matrix-vector products
        if 'gcon_prod' in evaluate or 'gcon_T_prod' in evaluate:
            if numpy.linalg.norm(x-self.cache['x']) > eps:
                # Previous evaluated point is *different* than the
                # point requested for the derivative. Recusively call
                # the routine with ['fobj', and 'fcon']
                self._masterFunc2(x, ['fobj', 'fcon'], writeHist=False)

            # Now, the point has been evaluated correctly so we
            # determine if we have to run the sens calc:
            if self.cache['gcon'] is None:
                timeA = time.time()
                gobj, gcon, fail = self.sens(
                    xuser, self.cache['fobj'], self.cache['fcon_user'])
                self.userSensTime += time.time()-timeA
                self.userSensCalls += 1
                # User values are stored immediately
                self.cache['gobj_user'] = copy.deepcopy(gobj)
                self.cache['gcon_user'] = copy.deepcopy(gcon)

                # Process objective gradient for optimizer
                gobj = self.optProb.processObjectiveGradient(gobj)

                # Process constraint gradients for optimizer
                gcon = self.optProb.processConstraintJacobian(gcon)
                gcon = self._convertJacobian(gcon)

                # Set cache values
                self.cache['gobj'] = gobj.copy()
                self.cache['gcon'] = gcon.copy()

                # Update fail flag
                masterFail = masterFail or fail
            else:
                # gcon is now in the cache, so we just retrieve it
                gcon = self.cache['gcon']

            if 'gcon_prod' in evaluate:
                outvec = gcon.dot(invec)
            else:
                outvec = gcon.T.dot(invec)

            # Set return values and history logging
            returns.append(outvec)
            hist['invec'] = invec
            hist['outvec'] = outvec

        elif 'jac_prod' in evaluate or 'jac_T_prod' in evaluate:
            # Call the true matrix-vector product functions, no matrix formation
            # Also, scale vectors appropriately before and after product
            timeA = time.time()
            if 'jac_prod' in evaluate:
                invec = self.optProb.invXScale.dot(invec)
                outvec, fail = self.sens[1](xuser, invec, sparse_only=sparse_only)
                outvec = self.optProb.conScale*outvec
            else:
                invec = self.optProb.conScale*invec
                outvec, fail = self.sens[2](xuser, invec, sparse_only=sparse_only)
                outvec = self.optProb.invXScale.dot(outvec)
            self.userSensTime += time.time() - timeA
            self.userSensCalls += 1

            # Update fail flag
            masterFail = masterFail or fail

            # Set return values and history logging
            returns.append(outvec)
            hist['invec'] = invec
            hist['outvec'] = outvec

        # end if 

        # Put the fail flag in the history
        hist['fail'] = masterFail

        # Write history if necessary
        if self.optProb.comm.rank == 0 and writeHist and self.storeHistory:
            self.hist.write(self.callCounter, hist)

        # We can now safely increment the call counter
        self.callCounter += 1

        # Tack the fail flag on at the end
        returns.append(masterFail)

        return returns


    def _on_setOption(self, name, value):
        """
        Set Optimizer Option Value (Optimizer Specific Routine)

        Parameters
        ----------
        name -> STRING: Option Name
        value: Option value
        """

        self.set_options.append([name, value])


    def _on_getOption(self, name):
        """
        Get Optimizer Option Value (Optimizer Specific Routine)

        Parameters
        ----------
        name -> STRING: Option Name
        """

        pass


    def _on_getInform(self, infocode):
        """
        Get Optimizer Result Information (Optimizer Specific Routine)

        Parameters
        ----------
        infocode -> INT: Status code
        """

        pass
