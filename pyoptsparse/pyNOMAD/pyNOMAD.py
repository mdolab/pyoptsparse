# File: pyNOMAD.py
# Author: Nicolas Bons
# Description: This program finds the optimum solution for a given objective
# function using the optimization package NOMAD.

from __future__ import absolute_import
from __future__ import print_function
# =============================================================================
# NOMAD Library
# =============================================================================
from . import nomad
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
# SLSQP Optimizer Class
# =============================================================================

class NOMAD(Optimizer):

    def __init__(self, *args, **kwargs):
        name = 'NOMAD'
        category = 'MADS Optimizer'
        defOpts = {
            # SLSQP Options
            'ACC': [float, 1e-6],         # Convergence Accurancy
            'MAXIT': [int, 1000000],          # Maximum Iterations
            'IPRINT': [int, 1],           # Output Level (<0 - None, 0 - Screen, 1 - File)
            'IOUT': [int, 6],             # Output Unit Number
            'IFILE': [str, 'NOMAD.out'],  # Output File Name
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
    def __call__(self, optProb, sens=None, sensStep=None, sensMode=None,
                 storeHistory=None, hotStart=None, storeSens=True):
        
        self.callCounter = 0
        self.storeSens = storeSens
        
        if len(optProb.constraints) == 0:
            # If the user *actually* has an unconstrained problem,
            # slsqp sort of chokes with that....it has to have at
            # least one constraint. So we will add one
            # automatically here:
            self.unconstrained = True
            optProb.dummyConstraint = False
        
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
            # Set history/hotstart
            self._setHistory(storeHistory, hotStart)
        
            # Define the objective function
            #--------------------------------------------------------------
            # x - list of design variables
            # fx - list of objective and constraints; first entry is objective

            def objfun(o, x_tuple):
                x = numpy.asarray(x_tuple)
                fobj, fcon, fail = self._masterFunc(x, ['fobj', 'fcon'])
                f = [fobj]
                g = fcon.tolist()
                funcallreturns = f + g
                return funcallreturns

            # Set objfun as the callback function
            problemset = nomad.NomadLinker()
            problemset.setCallback(objfun)

            # Setup argument list values
            
            acc = numpy.array([self.getOption('ACC')], numpy.float)
            maxit = self.getOption('MAXIT')
            iprint = self.getOption('IPRINT')
            ifile = self.getOption('IFILE')
            if iprint >= 0:
                if os.path.isfile(ifile):
                    os.remove(ifile)
            
            # Stopping criteria
            min_poll_size = 1e-6
            min_mesh_size = 0
                        
            # NOMAD Display (0-3)
            display_degree = 0      # 0 - no display
                        
            # Run NOMAD            
            t0 = time.time()
            solutionset = problemset.call(n, m+1, xs, blx, bux, min_poll_size, min_mesh_size, maxit, display_degree, iprint)
            ff = solutionset[0]
            optTime = time.time() -t0
            
            if self.storeHistory:
                self.hist.close()

            # Broadcast a -1 to indcate NOMAD has finished
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
