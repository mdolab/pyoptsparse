# File: pyNOMAD.py
# Author: Nicolas Bons
# Description: This program finds the optimum solution for a given objective
# function using the optimization package NOMAD.

from __future__ import absolute_import
from __future__ import print_function
# =============================================================================
# NOMAD Library
# =============================================================================
try:
    from . import nomad
except ImportError:
    nomad = None
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
            # NOMAD Options
            'maxiter':[int, 1000000],                   # Maximum number of function evaluations
            'minmeshsize':[float, 1e-12],           # Minimum refinement size of mesh
            'minpollsize':[float, 1e-12],           # Minimum step size for polling procedure
            'displaydegree':[int, 0],               # 0-none, 1-minimal display, 2-normal display, 3-full display
            'printfile':[int, 1]                    # 0-no output file, 1-output to NOMAD.out
            }
        informs = {}
        self.set_options = []
        if nomad is None:
            raise Error('There was an error importing the compiled'
                        'nomad module')

        Optimizer.__init__(self, name, category, defOpts, informs, *args, **kwargs)

    def __call__(self, optProb, sens=None, sensStep=None, sensMode=None,
                 storeHistory=None, hotStart=None, storeSens=True):

        if len(optProb.constraints) == 0:
            # If the user *actually* has an unconstrained problem,
            # slsqp sort of chokes with that....it has to have at
            # least one constraint. So we will add one
            # automatically here:
            self.unconstrained = True
            optProb.dummyConstraint = False

        self.callCounter = 0
        self.optProb = optProb
        self.optProb.finalizeDesignVariables()
        self.optProb.finalizeConstraints()
        self._setInitialCacheValues()
        blx, bux, xs = self._assembleContinuousVariables()
        xs = numpy.maximum(xs, blx)
        xs = numpy.minimum(xs, bux)
        n = len(xs)
        ff = self._assembleObjective()

        oneSided = True
        if self.unconstrained:
            m = 0
            meq = 0
        else:
            indices, blc, buc, fact = self.optProb.getOrdering(
                ['ne', 'le', 'ni', 'li'], oneSided=oneSided, noEquality=True)
            m = len(indices)

            self.optProb.jacIndices = indices
            self.optProb.fact = fact
            self.optProb.offset = buc

        if self.optProb.comm.rank == 0:
            # Set history/hotstart
            self._setHistory(storeHistory, hotStart)

            # Define the objective function
            #--------------------------------------------------------------
            def objfun(o, x_tuple):
                x = numpy.asarray(x_tuple)
                fail = False
                fobj, fcon, fail = self._masterFunc(x, ['fobj', 'fcon'])
                f = [fobj]
                g = fcon.tolist()
                if not fail:
                    cnt_eval = [0]
                else:
                    cnt_eval = [1]
                funcallreturns = cnt_eval + f + g
                return funcallreturns

            # Set objfun as the callback function
            problemset = nomad.NomadLinker()
            problemset.setCallback(objfun)

            # Setup NOMAD options
            maxit = self.getOption('maxiter')
            min_mesh_size = self.getOption('minmeshsize')
            min_poll_size = self.getOption('minpollsize')
            display_degree = self.getOption('displaydegree')
            print_file = self.getOption('printfile')

            # Run NOMAD
            t0 = time.time()
            solutionset = problemset.call(n, m+2, xs, blx, bux, min_poll_size, min_mesh_size, maxit, display_degree, print_file)
            ff = solutionset[1]
            optTime = time.time() -t0

            if self.storeHistory:
                self.hist.close()

            # Broadcast a -1 to indicate NOMAD has finished
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
