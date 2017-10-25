#/bin/env python
"""
pyNLPY_AUGLAG - A wrapper for the matrix-free augmented Lagrangian algorithm 
from the NLPy package. This optimizer can exploit sparsity in the constraint 
blocks.

Some assumptions to get this wrapper working:
- the NLPy package has been downloaded from Github and installed locally

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
try:
    from nlpy.model.mfnlp import MFModel
    from nlpy.optimize.solvers.sbmin import SBMINTotalLqnFramework
    from nlpy.optimize.solvers.sbmin import SBMINPartialLqnFramework
    from nlpy.optimize.solvers.sbmin import SBMINSplitLqnFramework
    from nlpy.optimize.solvers.tron import TronPartialLqnFramework
    from nlpy.optimize.solvers.tron import TronSplitLqnFramework
    from nlpy.optimize.solvers.tron import TronTotalLqnFramework
    from nlpy.optimize.solvers.auglag2 import AugmentedLagrangianPartialLsr1Framework
    from nlpy.optimize.solvers.auglag2 import AugmentedLagrangianPartialLsr1TronFramework
    from nlpy.optimize.solvers.auglag2 import AugmentedLagrangianSplitLsr1Framework
    from nlpy.optimize.solvers.auglag2 import AugmentedLagrangianSplitSr1Framework
    from nlpy.optimize.solvers.auglag2 import AugmentedLagrangianSplitLsr1TronFramework
    from nlpy.optimize.solvers.auglag2 import AugmentedLagrangianTotalLsr1AdjBroyAFramework
    from nlpy.optimize.solvers.auglag2 import AugmentedLagrangianTotalSr1AdjBroyAFramework
    from nlpy.optimize.solvers.auglag2 import AugmentedLagrangianTotalLsr1AdjBroyATronFramework
except:
    MFModel=None
# =============================================================================
# Standard Python modules
# =============================================================================
import os
import time
import copy
import signal
from contextlib import contextmanager
# =============================================================================
# External Python modules
# =============================================================================
import numpy
eps = numpy.finfo(1.0).eps
import logging
# ===========================================================================
# Extension modules
# ===========================================================================
from ..pyOpt_optimizer import Optimizer
from ..pyOpt_gradient import Gradient
from ..pyOpt_error import Error

# =============================================================================
# Timeout Class
# =============================================================================
class TimeoutException(Exception): pass

@contextmanager
def time_limit(seconds):
    def signal_handler(signum, frame):
        raise TimeoutException("Timed out!")
    signal.signal(signal.SIGALRM, signal_handler)
    signal.alarm(seconds)
    try:
        yield
    finally:
        signal.alarm(0)

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
        'Save Current Point':[bool,True],
        'Warm Restart':[bool,False],
        'Absolute Optimality Tolerance':[float,1.0e-6],
        'Relative Optimality Tolerance':[float,1.0e-6],
        'Absolute Feasibility Tolerance':[float,1.0e-6],
        'Relative Feasibility Tolerance':[float,1.0e-6],
        'Number of Quasi-Newton Pairs':[int,5],
        'Feasibility Quasi-Newton Pairs':[int,5],
        'Use N-Y Backtracking':[bool,True],
        'Use Magical Steps':[bool,True],
        # 'Use Tron':[bool,False],
        'Use Least-Squares Multipliers':[bool,False],
        'Use Damped Multiplier Update':[bool,True],
        'Use Quasi-Newton Jacobian':[bool,True],
        'Use Limited-Memory Approach':[bool,False],
        'Use Full-Memory Approach':[bool,False],
        'Use Tron':[bool,True],
        'Penalty Parameter':[float,10.],
        'Penalty Scaling':[float,0.1],
        'Maximum Inner Iterations':[int, 500],
        'Maximum Outer Iterations':[int, 20],
        'Maximum Time':[int,172000]
        }
        # Inform/Status codes go here
        informs = {
        1: 'Current point could not be improved',
        0: 'Successfully converged',
        -1: 'Maximum number of iterations reached',
        -2: 'Problem appears to be infeasible',
        -3: 'Solver stopped on user request',
        -5: 'Maximum run time exceeded'
        }
        self.set_options = []
        Optimizer.__init__(self, name, category, defOpts, informs, *args, **kwargs)

        # Additional Timing and call counters used by the matrix-free interface
        self.userJProdTime = 0.0
        self.userJTProdTime = 0.0
        self.userJProdCalls = 0
        self.userJTProdCalls = 0

    def __call__(self, optProb, sens=None, sensStep=None, sensMode=None,
                  storeHistory=None, hotStart=None, timeLimit=None,
                  storeSens=True):
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
            finite differenes use 'FD'. 'sens' may also be 'CS'
            which will cause pyOptSpare to compute the derivatives
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
            IDENTICAL**. As soon as he requested evaluation point does
            not match the history, function and gradient evaluations
            revert back to normal evaluations.

        timeLimit : number
            Number of seconds to run the optimization before a
            terminate flag is given to the optimizer and a "clean"
            exit is performed.

        storeSens : bool
            Flag sepcifying if sensitivities are to be stored in hist.
            This is necessay for hot-starting only.
            """

        self.callCounter = 0
        self.storeSens = storeSens
        
        if MFModel is None:
            raise Error('There was an error importing nlpy. nlpy must \
be installed to use NLPY_AUGLAG.')

        self.callCounter = 0

        # NLPy *can* handle unconstrained problems cleanly
        # Check if we still need this later
        # if len(optProb.constraints) == 0:
        #     self.unconstrained = True
        #     optProb.dummyConstraint = True

        # Save the optimization problem and finalize constraint
        # jacobian, in general can only do on root proc
        self.optProb = optProb
        self.optProb.finalizeDesignVariables()
        self.optProb.finalizeConstraints()

        self._setInitialCacheValues()
        blx, bux, xs = self._assembleContinuousVariables()

        # Redefined _setSens function for the matrix-free case
        self._setSens(sens, sensStep, sensMode)
        ff = self._assembleObjective()

        if self.optProb.nCon > 0:
            # We need to reorder this full jacobian...so get ordering:
            indices, blc, buc, fact = self.optProb.getOrdering(
                ['ne','ni','le','li'], oneSided=False)
            self.optProb.jacIndices = indices
            # Construct the inverse operator to self.optProb.jacIndices
            self.optProb.jacIndicesInv = numpy.argsort(self.optProb.jacIndices)
            self.optProb.fact = fact
            self.optProb.offset = numpy.zeros_like(indices)
            ncon = len(indices)

            # Count the number of nonlinear constraints for the quasi-Newton block
            # (A somewhat costly way to do this)
            tmp0, tmp1, tmp2, tmp3 = self.optProb.getOrdering(
                ['ne','ni'], oneSided=False)
            sparse_index = len(tmp0)
        else:
            blc = numpy.array([])
            buc = numpy.array([])
            sparse_index = 0

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

            # Configure tron logger for the case of using tron
            tronlogger = logging.getLogger('nlpy.tron')
            tronlogger.setLevel(logging.DEBUG)
            tronlogger.addHandler(hndlr)
            tronlogger.addHandler(hndlr2)
            tronlogger.propagate = False

            # This logger is only used for debugging
            # sr1logger = logging.getLogger('nlpy.lsr1')
            # sr1logger.setLevel(logging.DEBUG)
            # sr1logger.addHandler(hndlr)
            # sr1logger.propagate = False

            # pyOpt History logging - only do this on one processor?
            # dummy = self._setHistory(storeHistory,hotStart,coldStart,None)
            self._setHistory(storeHistory, hotStart)

        else:
            # Only the root proc stores the history
            self._setHistory("", hotStart)
        # end if

        # This optimizer has no hot start capability (too many vectors)
        # self._hotStart(storeHistory, hotStart)

        # Step 3: Pass options and solve the problem
        timeA = time.time()

        # Also need to pass the number of dense constraints
        # Assume the dense block is listed first in the problem definition
        # ** Simple sol'n: assume dense block is all nonlinear constraints **
        if self.options['Use Tron'][1]:
            if self.options['Use Quasi-Newton Jacobian'][1]:
                solver = AugmentedLagrangianTotalLsr1AdjBroyATronFramework(nlpy_problem, 
                    TronTotalLqnFramework, 
                    omega_abs=self.options['Absolute Optimality Tolerance'][1], 
                    eta_abs=self.options['Absolute Feasibility Tolerance'][1], 
                    omega_rel=self.options['Relative Optimality Tolerance'][1],
                    eta_rel=self.options['Relative Feasibility Tolerance'][1],
                    qn_pairs=self.options['Number of Quasi-Newton Pairs'][1],
                    least_squares_pi=self.options['Use Least-Squares Multipliers'][1],
                    data_prefix=self.options['Prefix'][1],
                    save_data=self.options['Save Current Point'][1],
                    warmstart=self.options['Warm Restart'][1],
                    sparse_index=sparse_index,
                    rho_init=self.options['Penalty Parameter'][1],
                    tau=self.options['Penalty Scaling'][1],
                    max_inner_iter=self.options['Maximum Inner Iterations'][1],
                    max_outer_iter=self.options['Maximum Outer Iterations'][1],
                    damped_pi=self.options['Use Damped Multiplier Update'][1])
            elif self.options['Use Limited-Memory Approach'][1]:
                solver = AugmentedLagrangianSplitLsr1TronFramework(nlpy_problem, 
                    TronSplitLqnFramework, 
                    omega_abs=self.options['Absolute Optimality Tolerance'][1], 
                    eta_abs=self.options['Absolute Feasibility Tolerance'][1], 
                    omega_rel=self.options['Relative Optimality Tolerance'][1],
                    eta_rel=self.options['Relative Feasibility Tolerance'][1],
                    qn_pairs=self.options['Number of Quasi-Newton Pairs'][1],
                    least_squares_pi=self.options['Use Least-Squares Multipliers'][1],
                    feas_qn_pairs=self.options['Feasibility Quasi-Newton Pairs'][1],
                    data_prefix=self.options['Prefix'][1],
                    save_data=self.options['Save Current Point'][1],
                    warmstart=self.options['Warm Restart'][1],
                    sparse_index=sparse_index,
                    rho_init=self.options['Penalty Parameter'][1],
                    tau=self.options['Penalty Scaling'][1],
                    max_inner_iter=self.options['Maximum Inner Iterations'][1],
                    max_outer_iter=self.options['Maximum Outer Iterations'][1],
                    damped_pi=self.options['Use Damped Multiplier Update'][1])
            else:
                # Try matrix-vector products with the exact Jacobian
                # Useful for comparisons, but not recommended for larger problems
                solver = AugmentedLagrangianPartialLsr1TronFramework(nlpy_problem, 
                    TronPartialLqnFramework, 
                    omega_abs=self.options['Absolute Optimality Tolerance'][1], 
                    eta_abs=self.options['Absolute Feasibility Tolerance'][1], 
                    omega_rel=self.options['Relative Optimality Tolerance'][1],
                    eta_rel=self.options['Relative Feasibility Tolerance'][1],
                    qn_pairs=self.options['Number of Quasi-Newton Pairs'][1],
                    least_squares_pi=self.options['Use Least-Squares Multipliers'][1],
                    data_prefix=self.options['Prefix'][1],
                    save_data=self.options['Save Current Point'][1],
                    warmstart=self.options['Warm Restart'][1],
                    rho_init=self.options['Penalty Parameter'][1],
                    tau=self.options['Penalty Scaling'][1],
                    max_inner_iter=self.options['Maximum Inner Iterations'][1],
                    max_outer_iter=self.options['Maximum Outer Iterations'][1],
                    damped_pi=self.options['Use Damped Multiplier Update'][1])
        else:
            if self.options['Use Quasi-Newton Jacobian'][1]:
                solver = AugmentedLagrangianTotalLsr1AdjBroyAFramework(nlpy_problem, 
                    SBMINTotalLqnFramework, 
                    omega_abs=self.options['Absolute Optimality Tolerance'][1], 
                    eta_abs=self.options['Absolute Feasibility Tolerance'][1], 
                    omega_rel=self.options['Relative Optimality Tolerance'][1],
                    eta_rel=self.options['Relative Feasibility Tolerance'][1],
                    qn_pairs=self.options['Number of Quasi-Newton Pairs'][1],
                    least_squares_pi=self.options['Use Least-Squares Multipliers'][1],
                    data_prefix=self.options['Prefix'][1],
                    save_data=self.options['Save Current Point'][1],
                    warmstart=self.options['Warm Restart'][1],
                    sparse_index=sparse_index,
                    rho_init=self.options['Penalty Parameter'][1],
                    tau=self.options['Penalty Scaling'][1],
                    max_inner_iter=self.options['Maximum Inner Iterations'][1],
                    max_outer_iter=self.options['Maximum Outer Iterations'][1],
                    damped_pi=self.options['Use Damped Multiplier Update'][1])
            elif self.options['Use Limited-Memory Approach'][1]:
                solver = AugmentedLagrangianSplitLsr1Framework(nlpy_problem, 
                    SBMINSplitLqnFramework, 
                    omega_abs=self.options['Absolute Optimality Tolerance'][1], 
                    eta_abs=self.options['Absolute Feasibility Tolerance'][1], 
                    omega_rel=self.options['Relative Optimality Tolerance'][1],
                    eta_rel=self.options['Relative Feasibility Tolerance'][1],
                    qn_pairs=self.options['Number of Quasi-Newton Pairs'][1],
                    least_squares_pi=self.options['Use Least-Squares Multipliers'][1],
                    feas_qn_pairs=self.options['Feasibility Quasi-Newton Pairs'][1],
                    data_prefix=self.options['Prefix'][1],
                    save_data=self.options['Save Current Point'][1],
                    warmstart=self.options['Warm Restart'][1],
                    sparse_index=sparse_index,
                    rho_init=self.options['Penalty Parameter'][1],
                    tau=self.options['Penalty Scaling'][1],
                    max_inner_iter=self.options['Maximum Inner Iterations'][1],
                    max_outer_iter=self.options['Maximum Outer Iterations'][1],
                    damped_pi=self.options['Use Damped Multiplier Update'][1])
            elif self.options['Use Full-Memory Approach'][1]:
                solver = AugmentedLagrangianSplitSr1Framework(nlpy_problem,
                    SBMINSplitLqnFramework, 
                    omega_abs=self.options['Absolute Optimality Tolerance'][1], 
                    eta_abs=self.options['Absolute Feasibility Tolerance'][1], 
                    omega_rel=self.options['Relative Optimality Tolerance'][1],
                    eta_rel=self.options['Relative Feasibility Tolerance'][1],
                    qn_pairs=self.options['Number of Quasi-Newton Pairs'][1],
                    least_squares_pi=self.options['Use Least-Squares Multipliers'][1],
                    feas_qn_pairs=self.options['Feasibility Quasi-Newton Pairs'][1],
                    data_prefix=self.options['Prefix'][1],
                    save_data=self.options['Save Current Point'][1],
                    warmstart=self.options['Warm Restart'][1],
                    sparse_index=sparse_index,
                    rho_init=self.options['Penalty Parameter'][1],
                    tau=self.options['Penalty Scaling'][1],
                    max_inner_iter=self.options['Maximum Inner Iterations'][1],
                    max_outer_iter=self.options['Maximum Outer Iterations'][1],
                    damped_pi=self.options['Use Damped Multiplier Update'][1])
            else:
                # Try matrix-vector products with the exact Jacobian
                # Useful for comparisons, but not recommended for larger problems
                solver = AugmentedLagrangianPartialLsr1Framework(nlpy_problem, 
                    SBMINPartialLqnFramework, 
                    omega_abs=self.options['Absolute Optimality Tolerance'][1], 
                    eta_abs=self.options['Absolute Feasibility Tolerance'][1], 
                    omega_rel=self.options['Relative Optimality Tolerance'][1],
                    eta_rel=self.options['Relative Feasibility Tolerance'][1],
                    qn_pairs=self.options['Number of Quasi-Newton Pairs'][1],
                    least_squares_pi=self.options['Use Least-Squares Multipliers'][1],
                    data_prefix=self.options['Prefix'][1],
                    save_data=self.options['Save Current Point'][1],
                    warmstart=self.options['Warm Restart'][1],
                    rho_init=self.options['Penalty Parameter'][1],
                    tau=self.options['Penalty Scaling'][1],
                    max_inner_iter=self.options['Maximum Inner Iterations'][1],
                    max_outer_iter=self.options['Maximum Outer Iterations'][1],
                    damped_pi=self.options['Use Damped Multiplier Update'][1])

        # if self.optProb.comm.rank == 0:
        #     print("Starting solve")
        try:
            with time_limit(self.options['Maximum Time'][1]):
                solver.solve(ny=self.options['Use N-Y Backtracking'], magic_steps_agg=self.options['Use Magical Steps'])
        except TimeoutException:
            solver.status = -5
            solver.tsolve = float(self.options['Maximum Time'][1])

        # Step 4: Collect and return solution
        optTime = time.time() - timeA

        if self.storeHistory:
            self.hist.close()

        sol_inform = {}
        sol_inform['value'] = solver.status
        sol_inform['text'] = self.informs[solver.status]
        xopt = solver.x[:self.optProb.ndvs].copy()
        sol = self._createSolution(optTime, sol_inform, solver.f, xopt)

        return sol


    def _clearTimings(self):
        """Clear timings and call counters"""
        Optimizer._clearTimings(self)
        self.userJProdTime = 0.0
        self.userJTProdTime = 0.0
        self.userJProdCalls = 0
        self.userJTProdCalls = 0


    def _setSens(self, sens, sensStep, sensMode):
        """
        For the matrix-free approach, the sens argument is actually a list 
        of three separate functions. The order of these functions must be 
        [obj_grad, jac_prod, jac_t_prod].

        Otherwise, this function is identical to that of the base class.
        """

        self.matrix_free = False
        if sens is None:
            raise Error("'None' value given for sens. Must be one \
of 'FD' or 'CS' or a user supplied function or group of functions.")
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
            raise Error("Unknown value given for sens. Must be None, 'FD', \
            'CS', a python function handle, or a list of handles")


    def _masterFunc3(self, x, invec, evaluate, writeHist=True, sparse_only=False):
        """
        A shell function for the matrix-free case, called on all processors.

        ** Right now, we assume that the 'evaluate' list has only one element **
        """

        timeAA = time.time()

        xscaled = self.optProb.invXScale*x
        xuser = self.optProb.processX(xscaled)

        masterFail = False

        # History storage is a little different for the matrix-free case.

        # Storing the whole set of input and output vectors is ridiculously 
        # expensive, so we only store function information.

        # Set basic parameters in history
        hist = {'xuser': xuser}
        returns = []

        # Evaluate the gradient of the objective function only
        if 'gobj' in evaluate:
            # if self.optProb.comm.rank == 0:
            # print("userObjCalls = %d"%self.userObjCalls)
            # print("userSensCalls = %d"%self.userSensCalls)
            if numpy.linalg.norm(x-self.cache['x']) > eps:
                # Previous evaluated point is *different* than the
                # point requested for the derivative. Recusively call
                # the routine with ['fobj', and 'fcon']
                self._masterFunc2(x, ['fobj', 'fcon'], writeHist=False)

            if self.cache['gobj'] is None:
                # Only the objective function is cached in this version
                timeA = time.time()
                gobj, fail = self.sens[0](xuser, funcs=self.cache['funcs'])
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
            # hist['gobj'] = self.cache['gobj_user']

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
                # gobj, gcon, fail = self.sens(
                #     xuser, self.cache['fobj'], self.cache['fcon_user'])
                funcsSens, fail = self.sens(xuser, self.cache['funcs'])
                self.userSensTime += time.time()-timeA
                self.userSensCalls += 1
                # User values are stored immediately
                # self.cache['gobj_user'] = copy.deepcopy(gobj)
                # self.cache['gcon_user'] = copy.deepcopy(gcon)
                self.cache['funcsSens'] = copy.deepcopy(funcsSens)

                # Process objective gradient for optimizer
                gobj = self.optProb.processObjectiveGradient(funcsSens)

                # Process constraint gradients for optimizer
                gcon = self.optProb.processConstraintJacobian(funcsSens)
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
            # hist['invec'] = invec
            # hist['outvec'] = outvec

        elif 'jac_prod' in evaluate or 'jac_T_prod' in evaluate:
            # Call the true matrix-vector product functions, no matrix formation
            # Convert input arrays to dictionaries before calling the callback function
            # Also, scale vectors appropriately before and after the products are formed
            if 'jac_prod' in evaluate:
                invec = self.optProb.invXScale*invec
                invec = self.optProb.processX(invec)
                timeA = time.time()
                outvec, fail = self.sens[1](xuser, invec, sparse_only=sparse_only, funcs=self.cache['funcs'])
                self.userJProdTime += time.time() - timeA
                self.userJProdCalls += 1
                outvec = self.optProb.processConstraints(outvec)
            else:
                invec = self.optProb.deProcessConstraints(invec)
                timeA = time.time()
                outvec, fail = self.sens[2](xuser, invec, sparse_only=sparse_only, funcs=self.cache['funcs'])
                self.userJTProdTime += time.time() - timeA
                self.userJTProdCalls += 1
                outvec = self.optProb.deProcessX(outvec)
                outvec = self.optProb.invXScale*outvec
                
            # prodTime = time.time() - timeA
            # if 'jac_prod' in evaluate:
            # else:

            # Update fail flag
            masterFail = masterFail or fail

            # Set return values and history logging
            returns.append(outvec)
            # hist['invec'] = invec
            # hist['outvec'] = outvec

        # end if 

        # Put the fail flag in the history
        # hist['fail'] = masterFail

        # Write history if necessary
        # if self.optProb.comm.rank == 0 and writeHist and self.storeHistory:
        #     self.hist.write(self.callCounter, hist)

        # We can now safely increment the call counter
        self.callCounter += 1

        # Tack the fail flag on at the end
        returns.append(masterFail)
        # Interface time specific to the Jacobian vector products
        self.interfaceTime += time.time() - timeAA

        return returns


    def _createSolution(self, optTime, sol_inform, obj, xopt):
        """
        Create the solution for the optimizer and append the data that is 
        specific to the matrix-free optimizer.
        """

        sol = Optimizer._createSolution(self, optTime, sol_inform, obj, xopt)
        sol.userJProdTime = self.userJProdTime
        sol.userJProdCalls = self.userJProdCalls
        sol.userJTProdTime = self.userJTProdTime
        sol.userJTProdCalls = self.userJTProdCalls

        # Recompute the interface time and optimizer code time to account for the 
        # separate matrix-vector product interface
        sol.interfaceTime = sol.interfaceTime - self.userJProdTime - self.userJTProdTime
        sol.optCodeTime = sol.optTime - self.interfaceTime

        # Since NLPy optimizers exploit MPI, we have to add the objective 
        # function here since the broadcasting of the base class is not done
        sol.objFun = self.optProb.objFun

        return sol


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
