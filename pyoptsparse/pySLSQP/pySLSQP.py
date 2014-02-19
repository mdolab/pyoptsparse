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
    from . import slsqp
except:
    raise ImportError('SLSQP shared library failed to import')
# =============================================================================
# Standard Python modules
# =============================================================================
import os
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
# =============================================================================
# SLSQP Optimizer Class
# =============================================================================
class SLSQP(Optimizer):
    """
    SLSQP Optimizer Class - Inherited from Optimizer Abstract Class
    """
    def __init__(self, *args, **kwargs):
        name = 'SLSQP'
        category = 'Local Optimizer'
        defOpts = {
            # SLSQP Options
            'ACC': [float, 1e-6],         # Convergence Accurancy
            'MAXIT': [int, 500],          # Maximum Iterations
            'IPRINT': [int, 1],           # Output Level (<0 - None, 0 - Screen, 1 - File)
            'IOUT': [int, 6],             # Output Unit Number
            'IFILE': [str, 'SLSQP.out'],  # Output File Name
            }
        informs = {
            -1 : "Gradient evaluation required (g & a)",
             0 : "Optimization terminated successfully.",
             1 : "Function evaluation required (f & c)",
             2 : "More equality constraints than independent variables",
             3 : "More than 3*n iterations in LSQ subproblem",
             4 : "Inequality constraints incompatible",
             5 : "Singular matrix E in LSQ subproblem",
             6 : "Singular matrix C in LSQ subproblem",
             7 : "Rank-deficient equality constraint subproblem HFTI",
             8 : "Positive directional derivative for linesearch",
             9 : "Iteration limit exceeded",
             }

        self.set_options = []
        Optimizer.__init__(self, name, category, defOpts, informs, *args,
                           **kwargs)

        # SLSQP needs jacobians in dense format
        self.jacType = 'dense2d'

    def __call__(self, optProb, sens=None, sensStep=None, sensMode=None,
                 storeHistory=None, hotStart=None, coldStart=None, 
                 timeLimit=None):
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
            from SLSQP does not match the history, function and
            gradient evaluations revert back to normal evaluations.

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
        # Set the number of nonlinear constraints snopt *thinks* we have:
        if self.unconstrained:
            m = 0
            meq = 0
        else:
            indices, blc, buc, fact = self.optProb.getOrdering(
                ['ne','le','ni','li'], oneSided=oneSided)
            m = len(indices)

            self.optProb.jacIndices = indices
            self.optProb.fact = fact
            self.optProb.offset = buc

            # Also figure out the number of equality:
            tmp0, __, __, __ = self.optProb.getOrdering(
                ['ne','le'], oneSided=oneSided)
            meq = len(tmp0)

        if self.optProb.comm.rank == 0:
            # Set history
            self._setHistory(storeHistory)

            if coldStart is not None:
                res1 = self._coldStart(coldStart)
                if res1 is not None:
                    xs[0:n] = res1

            #=================================================================
            # SLSQP - Objective/Constraint Values Function
            #=================================================================
            def slfunc(m, me, la, n, f, g, x):
                fobj, fcon, fail = self.masterFunc(x, ['fobj', 'fcon'])
                f = fobj
                g[0:m] = -fcon
                return f, g

            #=================================================================
            # SLSQP - Objective/Constraint Gradients Function
            #=================================================================
            def slgrad(m, me, la, n, f, g, df, dg, x):
                gobj, gcon, fail = self.masterFunc(x, ['gobj', 'gcon'])
                df[0:n] = gobj.copy()
                dg[0:m,0:n] = -gcon.copy()
                return df, dg

            # Setup hot start if necessary
            self._hotStart(storeHistory, hotStart)

            # Setup argument list values
            la = max(m, 1)
            gg = numpy.zeros([la], numpy.float)
            n1 = numpy.array([n+1], numpy.int)
            df = numpy.zeros([n+1], numpy.float)
            dg = numpy.zeros([la, n+1], numpy.float)
            acc = numpy.array([self.getOption('ACC')], numpy.float)
            maxit = self.getOption('MAXIT')
            iprint = self.getOption('IPRINT')
            iout = self.getOption('IOUT')
            ifile = self.getOption('IFILE')
            if iprint >= 0:
                if os.path.isfile(ifile):
                    os.remove(ifile)

            mode = 0
            mineq = m - meq + 2*(n+1)
            lsq = (n+1)*((n+1)+1) + meq*((n+1)+1) + mineq*((n+1)+1)
            lsi = ((n+1)-meq+1)*(mineq+2) + 2*mineq
            lsei = ((n+1)+mineq)*((n+1)-meq) + 2*meq + (n+1)
            slsqpb = (n+1)*(n/2) + 2*m + 3*n + 3*(n+1) + 1
            lwM = lsq + lsi + lsei + slsqpb + n + m
            lw = numpy.array([lwM], numpy.int)
            w = numpy.zeros([lw], numpy.float)
            ljwM = max(mineq, (n+1)-meq)
            ljw = numpy.array([ljwM], numpy.int)
            jw = numpy.zeros([ljw], numpy.intc)
            nfunc = numpy.array([0], numpy.int)
            ngrad = numpy.array([0], numpy.int)

            # Run SLSQP
            t0 = time.time()
            slsqp.slsqp(m, meq, la, n, xs, blx, bux, ff, gg, df, dg, acc, maxit,
                        iprint, iout, ifile, mode, w, lw, jw, ljw, nfunc,
                        ngrad, slfunc, slgrad)
            optTime = time.time() - t0

            if iprint > 0:
                slsqp.closeunit(self.getOption('IOUT'))

            if MPI:
                # Broadcast a -1 to indcate SLSQP has finished
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
            self.waitLoop()
            sol = None

        # Communication solution and return
        sol = self._communicateSolution(sol)

        return sol

    def _on_setOption(self, name, value):
        pass

    def _on_getOption(self, name, value):
        pass

