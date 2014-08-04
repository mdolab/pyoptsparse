#/bin/env python
"""
pyFSQP - A variation of the pyFSQP wrapper specificially designed to
work with sparse optimization problems.

Copyright (c) 2013-2014 by Dr. Gaetan Kenway
All rights reserved.

Tested on:
---------
Linux with intel

Developers:
-----------
- Dr. Gaetan Kenway (GKK)
History
-------
    v. 0.1    - Initial Wrapper Creation
"""
from __future__ import absolute_import
from __future__ import print_function
# =============================================================================
# FSQP Library
# =============================================================================
try:
    from . import ffsqp
except ImportError:
    ffsqp = None
# =============================================================================
# Standard Python modules
# =============================================================================
import os
import time
# =============================================================================
# External Python modules
# =============================================================================
import numpy
# # ===========================================================================
# # Extension modules
# # ===========================================================================
from ..pyOpt_optimizer import Optimizer
from ..pyOpt_error import Error
# =============================================================================
# FSQP Optimizer Class
# =============================================================================
class FSQP(Optimizer):
    """
    FSQP Optimizer Class - Inherited from Optimizer Abstract Class
    """
    def __init__(self, *args, **kwargs):
        
        name = 'FSQP'
        category = 'Local Optimizer'
        defOpts = {
        'mode':[int, 100],         # FSQP Mode (See Manual)
        'iprint':[int, 2],         # Output Level (0 - None, 1 - Final, 2 - Major, 3 - Major Details)
        'miter':[int, 500],        # Maximum Number of Iterations
        'bigbnd':[float, 1e10],    # Plus Infinity Value
        'epstol':[float, 1e-8],    # Convergence Tolerance
        'epseqn':[float, 0],       # Equality Constraints Tolerance
        'iout':[int, 6],           # Output Unit Number
        'ifile':[str, 'FSQP.out'], # Output File Name
        }
        informs = {
        0 : 'Normal termination of execution',
        1 : 'User-provided initial guess is infeasible for linear constraints, unable to generate a point satisfying all these constraints',
        2 : 'User-provided initial guess is infeasible for nonlinear inequality constraints and linear constraints, unable to generate a point satisfying all these constraints',
        3 : 'The maximum number of iterations has been reached before a solution is obtained',
        4 : 'The line search fails to find a new iterate',
        5 : 'Failure of the QP solver in attempting to construct d0, a more robust QP solver may succeed',
        6 : 'Failure of the QP solver in attempting to construct d1, a more robust QP solver may succeed',
        7 : 'Input data are not consistent, check print out error messages',
        8 : 'Two consecutive iterates are numerically equivalent before a stopping criterion is satisfied',
        9 : 'One of the penalty parameters exceeded bigbnd, the algorithm is having trouble satisfying a nonlinear equality constraint',
        }
        if ffsqp is None:
            raise Error('There was an error importing the compiled \
                        ffsqp module')

        Optimizer.__init__(self, name, category, defOpts, informs, *args, **kwargs)
        
        # We need jacobians in dens2d formation
        self.jacType = 'dense2d'
        
    def __call__(self, optProb, sens=None, sensStep=None, sensMode='FD',
                 storeHistory=None, hotStart=None, storeSens=True):
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
            IDENTICAL**. As soon as he requested evaluation point
            from SNOPT does not match the history, function and
            gradient evaluations revert back to normal evaluations.

        storeSens : bool
            Flag sepcifying if sensitivities are to be stored in hist.
            This is necessay for hot-starting only.
            """

        self.callCounter = 0
        self.storeSens = storeSens

        if len(optProb.constraints) == 0:
            self.unconstrained = True
            optProb.dummyConstraint = False

        self.optProb = optProb
        self.optProb.finalizeDesignVariables()
        self.optProb.finalizeConstraints()
              
        self._setInitialCacheValues()
        self._setSens(sens, sensStep, sensMode)
        blx, bux, xs = self._assembleContinuousVariables()
        ff = self._assembleObjective()
   
        # Determine all the constraint information, numbers etc. 
        if self.optProb.nCon > 0:
            # We need to reorder this full jacobian...so get ordering:
            indices, blc, buc, fact = self.optProb.getOrdering(
                ['ni','li','ne','le'], oneSided=True)
            ncon = len(indices)

            self.optProb.jacIndices = indices
            self.optProb.fact = fact
            self.optProb.offset = buc

            # We need to call getOrdering a few more times to get
            # the remaining sizes:
            indices, __, __, __ = self.optProb.getOrdering(
                ['ni'], oneSided=True)
            nineqn = len(indices)

            indices, __, __, __ = self.optProb.getOrdering(
                ['ni', 'li'], oneSided=True)
            nineq = len(indices)

            indices, __, __, __ = self.optProb.getOrdering(
                ['ne'], oneSided=True)
            neqn = len(indices)

            indices, __, __, __ = self.optProb.getOrdering(
                ['ne', 'le'], oneSided=True)
            neq = len(indices)
        else:
            nineqn = 0
            nineq = 0
            neqn = 0
            neq = 0
            ncon = 0
     
        # We make a split here: If the rank is zero we setup the
        # problem and run SNOPT, otherwise we go to the waiting loop:
        if self.optProb.comm.rank == 0:
            # Set history/hotstart/coldstart
            self._setHistory(storeHistory, hotStart)

            #======================================================================
            # FSQP - Objective Values Function
            #======================================================================
            def obj(nparam, j, x, fj):
                if self._checkEval(x):
                    self._internalEval(x)

                fj = self.storedData['fobj']

                return fj

            #======================================================================
            # FSQP - Constraint Values Function
            #======================================================================
            def cntr(nparam, j, x, gj):

                # for given j, assign to gj the value of the jth constraint evaluated at x
                if self._checkEval(x):
                    self._internalEval(x)

                gj = self.storedData['fcon'][j-1]

                return gj

            #======================================================================
            # FSQP - Objective Gradients Function
            #======================================================================
            def gradobj(nparam, j, x, gradfj, obj):

                # assign to gradfj the gradient of the jth objective function evaluated at x
                if self._checkEval(x):
                    self._internalEval(x)

                gradfj[0:nparam] = self.storedData['gobj']

                return gradfj

            #======================================================================
            # FSQP - Constraint Gradients Function
            #======================================================================
            def gradcntr(nparam, j, x, gradgj, obj):

                # assign to gradgj the gradient of the jth constraint evaluated at x
                if self._checkEval(x):
                    self._internalEval(x)

                gradgj[0:nparam] = self.storedData['gcon'][j-1]

                return gradgj

            # Setup argument list values
            nparam = len(xs)
            nvar = nparam
            nf = 1

            mode = self.getOption('mode')
            if self.getOption('iprint') >= 0:
                iprint = self.getOption('iprint')
            else:
                raise Error('Incorrect iprint option.  Must be >= 0')
            iout = self.getOption('iout')
            ifile = self.getOption('ifile')

            if iprint > 0:
                if os.path.isfile(ifile):
                    os.remove(ifile)

            gg = numpy.zeros(max(ncon,1))
            miter = self.getOption('miter')
            inform = 0
            bigbnd = self.getOption('bigbnd')
            epstol = self.getOption('epstol')
            epsneq = self.getOption('epseqn')
            udelta = 0
            nobj = 1
            iwsize = 6*nvar + 8*max([1,ncon]) + 7*max([1,nobj]) + 30
            iw = numpy.zeros([iwsize], numpy.float)
            nwsize = (4*nvar**2 + 5*max([1,ncon])*nvar + 3*max([1,nobj])*nvar + 
                      26*(nvar+max([1,nobj])) + 45*max([1,ncon]) + 100)
            w = numpy.zeros([nwsize], numpy.float)

            # Run FSQP
            t0 = time.time()
            ffsqp.ffsqp(nparam, nf, nineqn, nineq, neqn, neq, mode, iprint, miter,
                        inform, bigbnd, epstol, epsneq, udelta, blx, bux, xs, ff, 
                        gg, iw, iwsize, w, nwsize, obj, cntr, gradobj, gradcntr,
                        iout, ifile)
            optTime = time.time() - t0

            if iprint > 0:
                ffsqp.closeunit(iprint)
                
            # Broadcast a -1 to indcate SLSQP has finished
            self.optProb.comm.bcast(-1, root=0)

            # Store Results
            sol_inform = {}
            # sol_inform['value'] = inform
            # sol_inform['text'] = self.informs[inform[0]]

            # Create the optimization solution
            sol = self._createSolution(optTime, sol_inform, ff, xs)

        else:  # We are not on the root process so go into waiting loop:
            self._waitLoop()
            sol = None

        # Communication solution and return
        sol = self._communicateSolution(sol)

        return sol

    def _on_setOption(self, name, value):
        pass
        
    def _on_getOption(self, name):
        pass
