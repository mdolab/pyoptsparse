# /bin/env python
"""
pyALPSO - A pyOptSparse interface to ALPSO
work with sparse optimization problems.

Copyright (c) 2013-2014 by Dr. Gaetan Kenway
All rights reserved.

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
# Standard Python modules
# =============================================================================
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
# ALPSO Optimizer Class
# =============================================================================
class ALPSO(Optimizer):
    """
    ALPSO Optimizer Class - Inherited from Optimizer Abstract Class

    *Keyword arguments:**

    - pll_type -> STR: ALPSO Parallel Implementation (None, SPM- Static, DPM- Dynamic, POA-Parallel Analysis), *Default* = None
    """

    def __init__(self, *args, **kwargs):

        from . import alpso
        self.alpso = alpso

        category = 'Global Optimizer'
        defOpts = {
            'SwarmSize': [int, 40],  # Number of Particles (Depends on Problem dimensions)
            'maxOuterIter': [int, 200],  # Maximum Number of Outer Loop Iterations (Major Iterations)
            'maxInnerIter': [int, 6],  # Maximum Number of Inner Loop Iterations (Minor Iterations)
            'minInnerIter': [int, 6],  # Minimum Number of Inner Loop Iterations (Dynamic Inner Iterations)
            'dynInnerIter': [int, 0],  # Dynamic Number of Inner Iterations Flag
            'stopCriteria': [int, 1],  # Stopping Criteria Flag (0 - maxIters, 1 - convergence)
            'stopIters': [int, 5],  # Consecutive Number of Iterations for which the Stopping Criteria must be Satisfied
            'etol': [float, 1e-3],  # Absolute Tolerance for Equality constraints
            'itol': [float, 1e-3],  # Absolute Tolerance for Inequality constraints
            # 'ltol':[float, 1e-2],            # Absolute Tolerance for Lagrange Multipliers
            'rtol': [float, 1e-2],  # Relative Tolerance for Lagrange Multipliers
            'atol': [float, 1e-2],  # Absolute Tolerance for Lagrange Function
            'dtol': [float, 1e-1],  # Relative Tolerance in Distance of All Particles to Terminate (GCPSO)
            'printOuterIters': [int, 0],  # Number of Iterations Before Print Outer Loop Information
            'printInnerIters': [int, 0],  # Number of Iterations Before Print Inner Loop Information
            'rinit': [float, 1.0],  # Initial Penalty Factor
            'xinit': [int, 0],  # Initial Position Flag (0 - no position, 1 - position given)
            'vinit': [float, 1.0],  # Initial Velocity of Particles in Normalized [-1, 1] Design Space
            'vmax': [float, 2.0],  # Maximum Velocity of Particles in Normalized [-1, 1] Design Space
            'c1': [float, 2.0],  # Cognitive Parameter
            'c2': [float, 1.0],  # Social Parameter
            'w1': [float, 0.99],  # Initial Inertia Weight
            'w2': [float, 0.55],  # Final Inertia Weight
            'ns': [int, 15],  # Number of Consecutive Successes in Finding New Best Position of Best Particle Before Search Radius will be Increased (GCPSO)
            'nf': [int, 5],  # Number of Consecutive Failures in Finding New Best Position of Best Particle Before Search Radius will be Increased (GCPSO)
            'dt': [float, 1.0],  # Time step
            'vcrazy': [float, 1e-4], # Craziness Velocity (Added to Particle Velocity After Updating the Penalty Factors and Langangian Multipliers)
            'fileout': [int, 1],  # Flag to Turn On Output to filename
            'filename': [str, 'ALPSO.out'], # We could probably remove fileout flag if filename or fileinstance is given
            'seed': [float, 0],  # Random Number Seed (0 - Auto-Seed based on time clock)
            'HoodSize': [int, 40],  # Number of Neighbours of Each Particle
            'HoodModel': [str, 'gbest'], # Neighbourhood Model (dl/slring - Double/Single Link Ring, wheel - Wheel, Spatial - based on spatial distance, sfrac - Spatial Fraction)
            'HoodSelf': [int, 1],  # Selfless Neighbourhood Model (0 - Include Particle i in NH i, 1 - Don't Include Particle i)
            'Scaling': [int, 1],   # Design Variables Scaling Flag (0 - no scaling, 1 - scaling between [-1, 1])
            'parallelType': [str, ''],  # Type of parallelization ('' or 'EXT')
        }
        informs = {}
        Optimizer.__init__(self, 'ALPSO', category, defOpts, informs, *args, **kwargs)

    def __call__(self, optProb, storeHistory=None, **kwargs):
        """
        This is the main routine used to solve the optimization
        problem.

        Parameters
        ----------
        optProb : Optimization or Solution class instance
            This is the complete description of the optimization problem
            to be solved by the optimizer

        storeHistory : str
            File name of the history file into which the history of
            this optimization will be stored

        Notes
        -----
        The kwargs are there such that the sens= argument can be
        supplied (but ignored here in alpso)
        """
        # ======================================================================
        # ALPSO - Objective/Constraint Values Function
        # ======================================================================
        def objconfunc(x):
            fobj, fcon, fail = self._masterFunc(x, ['fobj', 'fcon'])
            return fobj, fcon

        # Save the optimization problem and finalize constraint
        # jacobian, in general can only do on root proc
        self.optProb = optProb
        self.optProb.finalizeDesignVariables()
        self.optProb.finalizeConstraints()
        self._setInitialCacheValues()

        if len(optProb.constraints) == 0:
            self.unconstrained = True

        xl, xu, xs = self._assembleContinuousVariables()
        xs = numpy.maximum(xs, xl)
        xs = numpy.minimum(xs, xu)
        n = len(xs)
        ff = self._assembleObjective()
        types = [0] * len(xs)
        oneSided = True
        if self.unconstrained:
            m = 0
            me = 0
        else:
            indices, blc, buc, fact = self.optProb.getOrdering(
                ['ne', 'le', 'ni', 'li'], oneSided=oneSided, noEquality=False)
            m = len(indices)

            self.optProb.jacIndices = indices
            self.optProb.fact = fact
            self.optProb.offset = buc
            indices, __, __, __ = self.optProb.getOrdering(
                ['ne', 'le'], oneSided=oneSided, noEquality=False)
            me = len(indices)

        nobj = 1


        if self.optProb.comm.rank == 0:
            # Set history/hotstart/coldstart
            self._setHistory(storeHistory, None)

            # Setup argument list values
            opt = self.getOption

            dyniI = self.getOption('dynInnerIter')
            if dyniI == 0:
                self.setOption('minInnerIter', opt('maxInnerIter'))

            if not opt('stopCriteria') in [0, 1]:
                raise Error('Incorrect Stopping Criteria Setting')

            if opt('fileout') not in [0, 1, 2, 3]:
                raise Error('Incorrect fileout Setting')

            if opt('seed') == 0:
                self.setOption('seed', time.time())

            # As far as I can tell, there is no need for this bulk attribute.
            # ALPSO calls the objconfunc iteratively for each particle in the
            # swarm, so we can deal with them one at a time, just as the other
            # optimizers.
            # self.optProb.bulk = opt('SwarmSize')

            # Run ALPSO
            t0 = time.time()
            opt_x, opt_f, opt_g, opt_lambda, nfevals, rseed = self.alpso.alpso(
                n, m, me, types, xs, xl, xu, opt('SwarmSize'), opt('HoodSize'),
                opt('HoodModel'), opt('maxOuterIter'), opt('maxInnerIter'),
                opt('minInnerIter'), opt('stopCriteria'), opt('stopIters'),
                opt('etol'), opt('itol'), opt('rtol'), opt('atol'), opt('dtol'),
                opt('printOuterIters'), opt('printInnerIters'), opt('rinit'),
                opt('vinit'), opt('vmax'), opt('c1'), opt('c2'), opt('w1'),
                opt('w2'), opt('ns'), opt('nf'), opt('vcrazy'), opt('fileout'),
                opt('filename'), None, None, opt('seed'),
                opt('Scaling'), opt('HoodSelf'), objconfunc)
            optTime = time.time() - t0

            # Broadcast a -1 to indcate NSGA2 has finished
            self.optProb.comm.bcast(-1, root=0)

            # Store Results
            sol_inform = {}
            # sol_inform['value'] = inform
            # sol_inform['text'] = self.informs[inform[0]]

            # Create the optimization solution
            sol = self._createSolution(optTime, sol_inform, opt_f, opt_x)
            for key in sol.objectives.keys():
                sol.objectives[key].value = opt_f
        else:  # We are not on the root process so go into waiting loop:
            self._waitLoop()
            sol = None

        # Communication solution and return
        sol = self._communicateSolution(sol)

        return sol

    def _on_setOption(self, name, value):
        if name == 'parallelType':
            value = value.upper()
            if value == 'EXT':
                try:
                    from . import alpso_ext
                    self.alpso = alpso_ext
                except ImportError:
                    raise ImportError('pyALPSO: ALPSO EXT shared library failed to import.')

            else:
                raise ValueError("parallel_type must be either '' or 'EXT'.")

    def _on_getOption(self, name, value):
        pass

    def _communicateSolution(self, sol):
        if sol is not None:
            sol.userObjCalls = self.optProb.comm.allreduce(sol.userObjCalls)
            sol.comm = None
        sol = self.optProb.comm.bcast(sol)
        sol.objFun = self.optProb.objFun
        sol.comm = self.optProb.comm

        return sol
