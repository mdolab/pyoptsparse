#/bin/env python
"""
pyCONMIN - A variation of the pyCONMIN wrapper specificially designed to
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
# CONMIN Library
# =============================================================================
try:
    from . import conmin
except ImportError:
    conmin = None
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
# =============================================================================
# CONMIN Optimizer Class
# =============================================================================
class CONMIN(Optimizer):
    """
    CONMIN Optimizer Class - Inherited from Optimizer Abstract Class
    """
    def __init__(self, *args, **kwargs):
        name = 'CONMIN'
        category = 'Local Optimizer'
        defOpts = {
            'ITMAX':[int, 1e4], # Maximum Number of Iterations
            'DELFUN':[float, 1e-6], # Objective Relative Tolerance
            'DABFUN':[float, 1e-6], # Objective Absolute Tolerance
            'ITRM':[int, 5],
            'NFEASCT':[int, 20],
            'IPRINT':[int, 4],  # Print Control (0 - None, 1 - Final, 2,3,4 - Debug)
            'IOUT':[int, 6], # Output Unit Number
            'IFILE':[str, 'CONMIN.out'], # Output File Name
        }
        informs = {}
        if conmin is None:
            raise Error('There was an error importing the compiled \
                        conmin module')

        self.set_options = []
        Optimizer.__init__(self, name, category, defOpts, informs, *args,
                           **kwargs)

        # CONMIN needs jacobians in dense format
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
            from CONMIN does not match the history, function and
            gradient evaluations revert back to normal evaluations.

        storeSens : bool
            Flag sepcifying if sensitivities are to be stored in hist.
            This is necessay for hot-starting only.
            """

        self.callCounter = 0
        self.storeSens = storeSens

        if len(optProb.constraints) == 0:
            # If the user *actually* has an unconstrained problem,
            # slsqp sort of chokes with that....it has to have at
            # least one constraint. So we will add one
            # automatically here:
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
        n = len(xs)
        ff = self._assembleObjective()

        oneSided = True
        noEquality = True
        if self.unconstrained:
            m = 0
        else:
            indices, blc, buc, fact = self.optProb.getOrdering(
                ['ne', 'le', 'ni', 'li'], oneSided=oneSided,
                noEquality=noEquality)
            m = len(indices)

            self.optProb.jacIndices = indices
            self.optProb.fact = fact
            self.optProb.offset = buc

        if self.optProb.comm.rank == 0:
            # Set history/hotstart/coldstart
            self._setHistory(storeHistory, hotStart)

            #=================================================================
            # CONMIN - Objective/Constraint Values Function
            #=================================================================
            def cnmnfun(n1, n2, x, f, g):
                fobj, fcon, fail = self._masterFunc(x[0:ndv], ['fobj', 'fcon'])
                f = fobj
                g[0:ncn] = fcon

                return f, g

            #=================================================================
            # CONMIN - Objective/Constraint Gradients Function
            #=================================================================
            def cnmngrad(n1, n2, x, f, g, ct, df, a, ic, nac):

                gobj, gcon, fail = self._masterFunc(x[0:ndv], ['gobj', 'gcon'])
                df[0:ndv] = gobj.copy()

                # Only assign the gradients for constraints that are
                # actually active:
                nac = 0
                for j in range(ncn):
                    if g[j] >= ct:
                        a[0:ndv, nac] = gcon[j, :]
                        ic[nac] = j + 1
                        nac += 1
                return df, a, ic, nac

            # Setup argument list values
            ndv = len(xs)
            ncn = m
            nn1 = ndv + 2
            nn2 = ncn + 2*ndv
            nn3 = max(nn2, ndv)
            nn4 = max(nn2, ndv)
            nn5 = 2*nn4
            gg = numpy.zeros(ncn, numpy.float)
            if self.getOption('IPRINT') >= 0 and self.getOption('IPRINT')  <= 4:
                iprint = self.getOption('IPRINT')
            else:
                raise Error('IPRINT option must be >= 0 and <= 4')

            iout = self.getOption('IOUT')
            ifile = self.getOption('IFILE')

            # Check if file exists and remove if necessary
            if iprint > 0:
                if os.path.isfile(ifile):
                    os.remove(ifile)

            itmax = self.getOption('ITMAX')
            delfun = self.getOption('DELFUN')

            #finit, ginit = cnmnfun([],[],xx,ff,gg)
            dabfun = self.getOption('DABFUN')

            itrm = self.getOption('ITRM')
            nfeasct = self.getOption('ITRM')
            nfdg = 1 # User will supply all gradients

            # Counters for functions and gradients
            nfun = 0
            ngrd = 0

            # Run CONMIN
            t0 = time.time()
            conmin.conmin(ndv, ncn, xs, blx, bux, ff, gg,
                          nn1, nn2, nn3, nn4, nn5,
                          iprint, iout, ifile, itmax, delfun, dabfun, itrm,
                          nfeasct, nfdg, nfun, ngrd, cnmnfun, cnmngrad)
            optTime = time.time() - t0

            if iprint > 0:
                conmin.closeunit(self.getOption('IOUT'))

            # Broadcast a -1 to indcate SLSQP has finished
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

