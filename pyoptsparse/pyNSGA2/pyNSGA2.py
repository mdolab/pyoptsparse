#/bin/env python
"""
pySLSQP - A variation of the pySLSQP wrapper specificially designed to
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
# SLSQP Library
# =============================================================================
try:
    from . import nsga2
except:
    nsga2 = None
# =============================================================================
# Standard Python modules
# =============================================================================
import os, copy
import time
# =============================================================================
# External Python modules
# =============================================================================
import numpy
from mpi4py import MPI
# ===========================================================================
# Extension modules
# ===========================================================================
from ..pyOpt_optimizer import Optimizer
from ..pyOpt_solution import Solution
from ..pyOpt_error import Error
# =============================================================================
# NSGA2 Optimizer Class
# =============================================================================
class NSGA2(Optimizer):
    """
    NSGA2 Optimizer Class - Inherited from Optimizer Abstract Class
    """
    def __init__(self, *args, **kwargs):

        name = 'NSGA-II'
        category = 'Global Optimizer'
        defOpts = {
        'PopSize': [int, 100],
        'maxGen': [int, 1000],
        'pCross_real': [float, 0.6],
        'pMut_real': [float, 0.2],
        'eta_c': [float, 10],
        'eta_m': [float, 20],
        'pCross_bin': [float, 0.0],
        'pMut_bin': [float, 0.0],
        'PrintOut': [int, 1],  # Flag to Turn On Output to filename (0 - , 1 - , 2 - )
        'seed': [float,0],      # Random Number Seed (0 - Auto-Seed based on time clock)
        'xinit': [int, 0],       # Use Initial Solution Flag (0 - random population, 1 - use given solution)
        }
        informs = {}
        Optimizer.__init__(self, name, category, defOpts, informs, *args, **kwargs)
        
        if nsga2 is None:
            raise Error('There was an error importing the compiled \
                        nsga2 module')

    def __call__(self, optProb, sens=None, sensStep=None, sensMode=None,
                 storeHistory=None, hotStart=None, coldStart=None):
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
            derivatives with finite differenes use \'FD\'. \'sens\'
            may also be \'CS\' which will cause pyOptSpare to compute
            the derivatives using the complex step method. Finally,
            \'sens\' may be a python function handle which is expected
            to compute the sensitivities directly. For expensive
            function evaluations and/or problems with large numbers of
            design variables this is the preferred method.

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
            IDENTICAL**. As soon as he requested evaluation point
            from NSGA2 does not match the history, function and
            gradient evaluations revert back to normal evaluations.

        coldStart : str
            Filename of the history file to use for "cold"
            restart. Here, the only requirment is that the number of
            design variables (and their order) are the same. Use this
            method if any of the optimization parameters have changed.
            """

        #======================================================================
        # NSGA-II - Objective/Constraint Values Function
        #======================================================================
        def objconfunc(nreal, nobj, ncon, x, f, g):
            xx = numpy.array(x)
            fobj, fcon, fail = self._masterFunc(xx, ['fobj','fcon'])

            f[0] = fobj
            g[0:ncon] = -fcon[0:ncon]
           
            return f, g

        self.callCounter = 0

        if len(optProb.constraints) == 0:
            # If the user *actually* has an unconstrained problem,
            # slsqp sort of chokes with that....it has to have at
            # least one constraint. So we will add one
            # automatically here:
            self.unconstrained = True
            optProb.dummyConstraint = False

        # Save the optimization problem and finialize constraint
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

        # Variables Handling
        n = len(xs)
        x = nsga2.new_doubleArray(n)
        xl = nsga2.new_doubleArray(n)
        xu = nsga2.new_doubleArray(n)
        for i in range(n):
            nsga2.doubleArray_setitem(x, i, xs[i])
            nsga2.doubleArray_setitem(xl, i, blx[i])
            nsga2.doubleArray_setitem(xu, i, bux[i])

        # Set the number of nonlinear constraints snopt *thinks* we have:
        if self.unconstrained:
            m = 0
        else:
            indices, blc, buc, fact = self.optProb.getOrdering(
                ['ne','le','ni','li'], oneSided=oneSided, noEquality=True)
            m = len(indices)

            self.optProb.jacIndices = indices
            self.optProb.fact = fact
            self.optProb.offset = buc

        g = nsga2.new_doubleArray(m)
        l = 1
        f = nsga2.new_doubleArray(1)

        if self.optProb.comm.rank == 0:
            # Set history
            self._setHistory(storeHistory)

            if coldStart is not None:
                res1 = self._coldStart(coldStart)
                if res1 is not None:
                    xs[0:n] = res1

            # Setup hot start if necessary
            self._hotStart(storeHistory, hotStart)

            # Setup argument list values
            nfeval = 0
            opt = self.getOption

            if self.getOption('PrintOut') >=0 and self.getOption('PrintOut') <=2:
                printout = self.getOption('PrintOut')
            else:
                raise Error('Incorrect option PrintOut')

            seed = self.getOption('seed')
            if seed == 0:
                seed = time.time()
        
            # Run NSGA-II
            nsga2.set_pyfunc(objconfunc)
            t0 = time.time()
            nsga2.nsga2(n, m, l, f, x, g, nfeval, xl, xu, opt('PopSize'), opt('maxGen'), 
                        opt('pCross_real'), opt('pMut_real'), opt('eta_c'), opt('eta_m'), 
                        opt('pCross_bin'), opt('pMut_bin'), printout,seed, opt('xinit'))
            optTime = time.time() - t0

            if MPI:
                # Broadcast a -1 to indcate NSGA2 has finished
                self.optProb.comm.bcast(-1, root=0)

            # Store Results
            sol_inform = {}
            #sol_inform['value'] = inform
            #sol_inform['text'] = self.informs[inform[0]]

            # Create the optimization solution
            sol = self._createSolution(optTime, sol_inform, ff)

            # Now set the x-values:
            i = 0
            for dvSet in sol.variables.keys():
                for dvGroup in sol.variables[dvSet]:
                    for var in sol.variables[dvSet][dvGroup]:
                        var.value = xs[i]
                        i += 1

            sol.fStar = ff

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

        # fstar = [0.]*l
        # for i in xrange(l):
        #     fstar[i] = nsga2.doubleArray_getitem(f,i)
        #     i += 1
        # #end
        
        # xstar = [0.]*n
        # for i in xrange(n):
        #     xstar[i] = nsga2.doubleArray_getitem(x,i)
        # #end
        
        # inform = {}
        
        # return fstar, xstar, {'fevals':nfeval,'time':sol_time,'inform':inform}
        
