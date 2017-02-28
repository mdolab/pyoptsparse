#/bin/env python
"""pyNLPQLP - A pyOptSparse wrapper for Schittkowski's NLPQLP
optimizion algorithm.

Copyright (c) 2013-2015 by Dr. Gaetan Kenway
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
# NLPQL Library
# =============================================================================
try:
    from . import nlpqlp
except ImportError:
    nlpqlp = None
# =============================================================================
# Standard Python modules
# =============================================================================
import os
import time
# =============================================================================
# External Python modules
# =============================================================================
import numpy
# ===========================================================================
# Extension modules
# ===========================================================================
from ..pyOpt_optimizer import Optimizer
from ..pyOpt_error import Error
eps = numpy.finfo(1.0).eps
# =============================================================================
# NLPQL Optimizer Class
# =============================================================================
class NLPQLP(Optimizer):
    """
    NLPQL Optimizer Class - Inherited from Optimizer Abstract Class
    """
    def __init__(self, *args, **kwargs):
        name = 'NLPQLP'
        category = 'Local Optimizer'
        defOpts = {
            # NLPQL Options
            'accuracy': [float, 1e-6],    # Convergence Accurancy
            'accuracyQP': [float, 1e-14], # Convergence Accurancy for QP
            'stepMin': [float, 1e-6],     # Minimum step length
            'maxFun': [int, 20],          # Maximum Number of Function Calls During Line Search
            'maxIt': [int, 500],          # Maximum Number of Iterations
            'maxNM': [int, 1],            # Maximum stack size for non-monotone line search
            'rho': [float, 1.0],          # Factor scaling identify for IFAIL=2
            'iPrint': [int, 2],           # Output Level (0 - None, 1 - Final, 2 - Major, 3 - Major/Minor, 4 - Full)
            'mode': [int, 0],             # NLPQL Mode (0 - Normal Execution, 1 to 18 - See Manual)
            'iOut': [int, 6],             # Output Unit Number
            'lMerit': [bool, True],       # Merit Function Type (True - L2 Augmented Penalty, False - L1 Penalty)
            'lQl': [bool, False],         # QP Subproblem Solver (True - Quasi-Newton, False - Cholesky)
            'iFile': [str,'NLPQLP.out'],  # Output File Name
        }
        informs = {
            -2 : 'Compute gradient values w.r.t. the variables stored in' \
            ' first column of X, and store them in DF and DG.' \
            ' Only derivatives for active constraints ACTIVE(J)=.TRUE. need to be computed.',
            -1 : 'Compute objective fn and all constraint values subject' \
            'the variables found in the first L columns of X, and store them in F and G.',
            0 : 'The optimality conditions are satisfied.', 
            1 : ' The algorithm has been stopped after MAXIT iterations.',
            2 : ' The algorithm computed an uphill search direction.',
            3 : ' Underflow occurred when determining a new approximation matrix' \
            'for the Hessian of the Lagrangian.',
            4 : 'The line search could not be terminated successfully.', 
            5 : 'Length of a working array is too short.' \
            ' More detailed error information is obtained with IPRINT>0',
            6 : 'There are false dimensions, for example M>MMAX, N>=NMAX, or MNN2<>M+N+N+2.',
            7 : 'The search direction is close to zero, but the current iterate is still infeasible.',
            8 : 'The starting point violates a lower or upper bound.',
            9 : 'Wrong input parameter, i.e., MODE, LDL decomposition in D and C' \
            ' (in case of MODE=1), IPRINT, IOUT',
            10 : 'Internal inconsistency of the quadratic subproblem, division by zero.',
            100 : 'The solution of the quadratic programming subproblem has been' \
            ' terminated with an error message and IFAIL is set to IFQL+100,' \
            ' where IFQL denotes the index of an inconsistent constraint.',
        }
        if nlpqlp is None:
            raise Error("There was an error importing the compiled "
                        "nlpqlp module")

        Optimizer.__init__(self, name, category, defOpts, informs, *args, **kwargs)
        # NLPQLP needs jacobians in dense format
        self.jacType = 'dense2d'

    def __call__(self, optProb, sens=None, sensStep=None, sensMode=None,
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
            Specifiy method to compute sensitivities. To
            explictly use pyOptSparse gradient class to do the
            derivatives with finite differenes use 'FD'. 'sens'
            may also be 'CS' which will cause pyOptSpare to compute
            the derivatives using the complex step method. Finally,
            'sens' may be a python function handle which is expected
            to compute the sensitivities directly. For expensive
            function evaluations and/or problems with large numbers of
            design variables this is the preferred method.

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
            from NLPQL does not match the history, function and
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

        # Save the optimization problem and finalize constraint
        # jacobian, in general can only do on root proc
        self.optProb = optProb
        self.optProb.finalizeDesignVariables()
        self.optProb.finalizeConstraints()
        self._setInitialCacheValues()
        self._setSens(sens, sensStep, sensMode)
        blx, bux, xs = self._assembleContinuousVariables()
        xs = numpy.maximum(xs, blx)
        xs = numpy.minimum(xs, bux)
        nvar = len(xs)
        ff = self._assembleObjective()

        oneSided = True
        # Set the number of nonlinear constraints snopt *thinks* we have:
        if self.unconstrained:
            m = 0
            meq = 0
        else:
            indices, blc, buc, fact = self.optProb.getOrdering(
                ['ne', 'le', 'ni', 'li'], oneSided=oneSided)
            m = len(indices)

            self.optProb.jacIndices = indices
            self.optProb.fact = fact
            self.optProb.offset = buc

            # Also figure out the number of equality:
            tmp0, __, __, __ = self.optProb.getOrdering(
                ['ne', 'le'], oneSided=oneSided)
            meq = len(tmp0)

        if self.optProb.comm.rank == 0:
            # Set history/hotstart/coldstart
            self._setHistory(storeHistory, hotStart)

            #=================================================================
            # NLPQL - Objective/Constraint Values Function (Real Valued) 
            #=================================================================
            def nlfunc(m, me, mmax, n, f, g, x, active, fail):
                fobj, fcon, fail = self._masterFunc(x, ['fobj', 'fcon'])
                f = fobj
                g[0:m] = -fcon
                return f, g, fail

            #=================================================================
            # NLPQL - Objective/Constraint Gradients Function
            #=================================================================
            def nlgrad(m, me, mmax, n, f, g, df, dg, x, active, wa):
                gobj, gcon, fail = self._masterFunc(x, ['gobj', 'gcon'])
                df[0:n] = gobj.copy()
                dg[0:m, 0:n] = -gcon.copy()
                return df, dg

            # setup argument list values

            np = 1 # We only allow a single "processor" ie we are
                   # actually running NLPQL (no P)

            # Set som basic sizes
            m = m
            me = meq
            mmax = max(1, m)

            n = nvar
            nmax = max(2, n+2)
            mnn2 = m+n+n+2

            # xs, ff, and gg have to have an extra dimension
            # associated with them for the NP. We will do this
            # correctly even though np is hard-coded to 1. 
            xs = numpy.array([xs]).T
            f = numpy.array([ff])
            g = numpy.zeros((mmax, np))
                
            df = numpy.zeros(nmax)
            dg = numpy.zeros((mmax, nmax))
            u = numpy.zeros(mnn2)
            c = numpy.zeros((nmax, nmax))
            d = numpy.zeros(nmax)
            go = self.getOption
            if go('iPrint') < 0 or go('iPrint') > 4:
                raise Error('Incorrect iPrint option. Must be >=0 and <= 4')
            
            if not(go('mode') >= 0 and go('mode') <= 18):
                raise Error('Incorrect mode option. Must be >= 0 and <= 18.')

            if os.path.isfile(go('iFile')):
                os.remove(go('iFile'))
            ifail = 0
            # Run NLPQL
            t0 = time.time()

            nlpqlp.wrapper(np, m, me, mmax, n, nmax, mnn2, xs, f, g, df, dg, u,
                           blx, bux, c, d, go('accuracy'), go('accuracyQP'), 
                           go('stepMin'), go('maxFun'), go('maxIt'), go('maxNM'),
                           go('rho'), go('mode'), go('iPrint'), go('iOut'), 
                           go('iFile'), ifail, go('lMerit'), go('lQl'), 
                           nlfunc, nlgrad)

            optTime = time.time() - t0
               
            # Broadcast a -1 to indcate NLPQL has finished
            self.optProb.comm.bcast(-1, root=0)

            # Store Results
            sol_inform = {}
            #sol_inform['value'] = inform
            #sol_inform['text'] = self.informs[inform[0]]

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

    def _on_getOption(self, name, value):
        pass


