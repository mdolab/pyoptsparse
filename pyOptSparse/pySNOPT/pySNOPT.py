from __future__ import absolute_import
from __future__ import print_function
#/bin/env python
"""
pySNOPT - A variation of the pySNOPT wrapper specificially designed to
work with sparse optimization problems.

Copyright (c) 2013-2014 by Dr. Gaetan Kenway
All rights reserved.

Tested on:
---------
Linux with intel

Developers:
-----------
- Dr. Gaetan Kenway (GKK)
- Dr. Graeme Kennedy (GJK)
History
-------
    v. 0.1    - Initial Wrapper Creation 
"""
# =============================================================================
# SNOPT Library
# =============================================================================
try:
    from . import snopt
except:
    raise ImportError('SNOPT shared library failed to import')
# =============================================================================
# Standard Python modules
# =============================================================================
import os
import sys
import copy
import time
import types
# =============================================================================
# External Python modules
# =============================================================================
import numpy
import shelve
from scipy import sparse
# # ===========================================================================
# # Extension modules
# # ===========================================================================
from ..pyOpt_optimizer import Optimizer
from ..pyOpt_history import History
from ..pyOpt_solution import Solution
from ..pyOpt_error import Error
# =============================================================================
# Misc Definitions
# =============================================================================
inf = 1e20  # define a value for infinity

# Try to import mpi4py and determine rank
try: 
    from mpi4py import MPI
    rank = MPI.COMM_WORLD.rank
except:
    rank = 0
    MPI = None

# =============================================================================
# SNOPT Optimizer Class
# =============================================================================
class SNOPT(Optimizer):
    """
    SNOPT Optimizer Class - Inherited from Optimizer Abstract Class
    """
    
    def __init__(self, *args, **kwargs):
        """
        SNOPT Optimizer Class Initialization
        """
        
        name = 'SNOPT'
        category = 'Local Optimizer'
        defOpts = {
        # SNOPT Printing Options
        'outputDirectory':[str, './'],      # Directory to put files in.     
        'Major print level':[int,1],                     # Majors Print (1 - line major iteration log)
        'Minor print level':[int,1],                     # Minors Print (1 - line minor iteration log)
        'Print file':[str,'SNOPT_print.out'],            # Print File Name (specified by subroutine snInit)
        'iPrint':[int,18],                                # Print File Output Unit (override internally in snopt?)
        'Summary file':[str,'SNOPT_summary.out'],        # Summary File Name (specified by subroutine snInit)
        'iSumm':[int,19],                                # Summary File Output Unit (override internally in snopt?)
        'Print frequency':[int,100],                    # Minors Log Frequency on Print File
        'Summary frequency':[int,100],                    # Minors Log Frequency on Summary File
        'Solution':[str,'Yes'],                            # Print Solution on the Print File
        'Suppress options listing':[type(None),None],    # (options are normally listed)
        'System information':[str,'No'],                # Print System Information on the Print File
        # SNOPT Problem Specification Options
        'Problem Type':[str,'Minimize'],                # ('Maximize': alternative over Minimize, 'Feasible point': alternative over Minimize or Maximize)
        'Objective row':[int,1],                         # (has precedence over ObjRow (snOptA))
        'Infinite bound':[float,1.0e+20],                # Infinite Bound Value
        # SNOPT Convergence Tolerances Options
        'Major feasibility tolerance':[float,1.0e-6],    # Target Nonlinear Constraint Violation
        'Major optimality tolerance':[float,1.0e-6],     # Target Complementarity Gap
        'Minor feasibility tolerance':[float,1.0e-6],     # For Satisfying the QP Bounds
        # SNOPT Derivative Checking Options
        'Verify level':[int,0],                            # Gradients Check Flag
        # SNOPT Scaling Options
        'Scale option':[int,1],                            # Scaling (1 - linear constraints and variables)
        'Scale tolerance':[float,0.9],                    # Scaling Tolerance
        'Scale Print':[type(None),None],                # Default: scales are not printed
        # SNOPT Other Tolerances Options
        'Crash tolerance':[float,0.1],                    # 
        'Linesearch tolerance':[float,0.9],             # smaller for more accurate search
        'Pivot tolerance':[float,3.7e-11],                 # epsilon^(2/3)
        # SNOPT QP subproblems Options
        'QPSolver':[str,'Cholesky'],                     # Default: Cholesky
        'Crash option':[int,3],                         # (3 - first basis is essentially triangular)
        'Elastic mode':[str,'No'],                        # (start with elastic mode until necessary)
        'Elastic weight':[float,1.0e+4],                 # (used only during elastic mode)
        'Iterations limit':[int,10000],                 # (or 20*ncons if that is more)
        'Partial price':[int,1],                         # (10 for large LPs)
        'Start':[str,'Cold'],                             # has precedence over argument start, ('Warm': alternative to a cold start)
        # SNOPT SQP method Options
        'Major iterations limit':[int,1000],             # or ncons if that is more
        'Minor iterations limit':[int,500],             # or 3*ncons if that is more
        'Major step limit':[float,2.0],                    # 
        'Superbasics limit':[int,None],                 # (n1 + 1, n1 = number of nonlinear variables)
        'Derivative level':[int,3],                        # (NOT ALLOWED IN snOptA)
        'Derivative option':[int,1],                    # (ONLY FOR snOptA)
        'Derivative linesearch':[type(None),None],        #
        'Nonderivative linesearch':[type(None),None],    #
        'Function precision':[float,3.0e-13],             # epsilon^0.8 (almost full accuracy)
        'Difference interval':[float,5.5e-7],             # Function precision^(1/2)
        'Central difference interval':[float,6.7e-5],    # Function precision^(1/3)
        'New superbasics limit':[int,99],                # controls early termination of QPs
        'Objective row':[int,1],                        # row number of objective in F(x)
        'Penalty parameter':[float,0.0],                 # initial penalty parameter
        'Proximal point method':[int,1],                # (1 - satisfies linear constraints near x0)
        'Reduced Hessian dimension':[int,2000],            # (or Superbasics limit if that is less)
        'Violation limit':[float,10.0],                     # (unscaled constraint violation limit)
        'Unbounded step size':[float,1.0e+18],            # 
        'Unbounded objective':[float,1.0e+15],            # 
        # SNOPT Hessian approximation Options
        'Hessian full memory':[type(None),None],         # default if n1 <= 75
        'Hessian limited memory':[type(None),None],     # default if n1 > 75
        'Hessian frequency':[int,999999],                 # for full Hessian (never reset)
        'Hessian updates':[int,10],                         # for limited memory Hessian
        'Hessian flush':[int,999999],                     # no flushing
        # SNOPT Frequencies Options
        'Check frequency':[int,60],                         # test row residuals ||Ax - sk||
        'Expand frequency':[int,10000],                    # for anti-cycling procedure
        'Factorization frequency':[int,50],                # 100 for LPs
        'Save frequency':[int,100],                         # save basis map
        # SNOPT LUSOL Options
        'LU factor tolerance':[float,3.99],             # for NP (100.0 for LP)
        'LU update tolerance':[float,3.99],             # for NP ( 10.0 for LP)
        'LU singularity tolerance':[float,3.2e-11],        # 
        'LU partial pivoting':[type(None),None],         # default threshold pivoting strategy
        'LU rook pivoting':[type(None),None],             # threshold rook pivoting
        'LU complete pivoting':[type(None),None],         # threshold complete pivoting
        # SNOPT Basis files Options
        'Old basis file':[int,0],                         # input basis map
        'New basis file':[int,0],                         # output basis map
        'Backup basis file':[int,0],                     # output extra basis map
        'Insert file':[int,0],                             # input in industry format
        'Punch file':[int,0],                             # output Insert data
        'Load file':[int,0],                             # input names and values
        'Dump file':[int,0],                             # output Load data
        'Solution file':[int,0],                         # different from printed solution
        # SNOPT Partitions of cw, iw, rw Options
        'Total character workspace':[int,500],          # lencw: 500
        'Total integer workspace':[int,None],            # leniw: 500 + 100 * (m+n) 
        'Total real workspace':[int,None],                # lenrw: 500 + 200 * (m+n)
        'User character workspace':[int,500],            # 
        'User integer workspace':[int,500],                # 
        'User real workspace':[int,500],                # 
        #SNOPT Miscellaneous Options
        'Debug level':[int,0],                             # (0 - Normal, 1 - for developers)
        'Timing level':[int,3],                            # (3 - print cpu times)
        }
        informs = {
        0 : 'finished successfully',
        1 : 'optimality conditions satisfied',
        2 : 'feasible point found',
        3 : 'requested accuracy could not be achieved',
        4 : 'weak QP minimizer',
        10 : 'the problem appears to be infeasible',
        11 : 'infeasible linear constraints',
        12 : 'infeasible linear equalities',
        13 : 'nonlinear infeasibilities minimized',
        14 : 'infeasibilities minimized',
        15 : 'infeasible linear constraints in QP subproblem',
        20 : 'the problem appears to be unbounded',
        21 : 'unbounded objective',
        22 : 'constraint violation limit reached',
        30 : 'resource limit error',
        31 : 'iteration limit reached',
        32 : 'major iteration limit reached',
        33 : 'the superbasics limit is too small',
        40 : 'terminated after numerical difficulties',
        41 : 'current point cannot be improved',
        42 : 'singular basis',
        43 : 'cannot satisfy the general constraints',
        44 : 'ill-conditioned null-space basis',
        50 : 'error in the user-supplied functions',
        51 : 'incorrect objective  derivatives',
        52 : 'incorrect constraint derivatives',
        53 : 'the QP Hessian is indefinite',
        54 : 'incorrect second derivatives',
        55 : 'incorrect derivatives',
        60 : 'undefined user-supplied functions',
        61 : 'undefined function at the first feasible point',
        62 : 'undefined function at the initial point',
        63 : 'unable to proceed into undefined region',
        70 : 'user requested termination',
        71 : 'terminated during function evaluation',
        72 : 'terminated during constraint evaluation',
        73 : 'terminated during objective evaluation',
        74 : 'terminated from monitor routine',
        80 : 'insufficient storage allocated',
        81 : 'work arrays must have at least 500 elements',
        82 : 'not enough character storage',
        83 : 'not enough integer storage',
        84 : 'not enough real storage',
        90 : 'input arguments out of range',
        91 : 'invalid input argument',
        92 : 'basis file dimensions do not match this problem',
        93 : 'the QP Hessian is indefinite',
        100 : 'finished successfully',
        101 : 'SPECS file read',
        102 : 'Jacobian structure estimated',
        103 : 'MPS file read',
        104 : 'memory requirements estimated',
        105 : 'user-supplied derivatives appear to be correct',
        106 : 'no derivatives were checked',
        107 : 'some SPECS keywords were not recognized',
        110 : 'errors while processing MPS data',
        111 : 'no MPS file specified',
        112 : 'problem-size estimates too small',
        113 : 'fatal error in the MPS file',
        120 : 'errors while estimating Jacobian structure',
        121 : 'cannot find Jacobian structure at given point',
        130 : 'fatal errors while reading the SP',
        131 : 'no SPECS file (iSpecs le 0 or iSpecs gt 99)',
        132 : 'End-of-file while looking for a BEGIN',
        133 : 'End-of-file while reading SPECS file',
        134 : 'ENDRUN found before any valid SPECS',
        140 : 'system error',
        141 : 'wrong no of basic variables',
        142 : 'error in basis package',
        142 : 'Problem dimensions are too large'
        }
        self.set_options = []
        Optimizer.__init__(self, name, category, defOpts, informs, *args, **kwargs)

        # Snopt doesn't required appending linearConstraint 
        self.appendLinearConstraints = False

        # Snopt need jacobians in csc format
        self.jacType = 'csc'

        # Constrained until we know otherwise :-)
        self.unconstrained = False
    def __call__(self, optProb, sens=None, sensStep=None, sensMode=None,
                  storeHistory=None, hotStart=None, warmStart=None,
                  coldStart=None, timeLimit=None, comm=None):
        """
        This is the main routine used to solve the optimization
        problem.

        Parameters
        ----------
        optProb : Optimization or Solution class instance
            This is the complete description of the optimization problem
            to be solved by the optimizer

        sens : str or python Function.
            Specifiy method to compute sensitivities. The default is
            None which will use SNOPT\'s own finite differences which
            are vastly superiour to the pyOptSparse implementation. To
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
            from SNOPT does not match the history, function and
            gradient evaluations revert back to normal evaluations. 
             
        warmStart : str
            File name of the history file to use for "warm"
            restart. The only difference between a "warm" and "cold"
            restart is that the state of the variables is used when
            initializing snopt. This information is stored in the
            history file from snopt only. If a warm start cannot be
            preformed for any reason, a cold start is attempted
            instead.

        coldStart : str
            Filename of the history file to use for "cold"
            restart. Here, the only requirment is that the number of
            design variables (and their order) are the same. Use this
            method if any of the optimization parameters have changed.

        timeLimit : number
            Number of seconds to run the optimization before a
            terminate flag is given to the optimizer and a "clean"
            exit is performed.

        comm : MPI Intra communicator
            Specifiy a MPI comm to use. Default is None. If mpi4py is not
            available, the serial mode will still work. if mpi4py *is*
            available, comm defaluts to MPI.COMM_WORLD. 
            """
        
        self.callCounter = 0
  
        if len(optProb.constraints) == 0:
            # If the user *actually* has an unconstrained problem,
            # snopt sort of chokes with that....it has to have at
            # least one constraint. So we will add one
            # automatically here:
            self.unconstrained = True
            optProb.dummyConstraint = True
            
        # Save the optimization problem and finialize constraint
        # jacobian, in general can only do on root proc

        self.optProb = optProb
        self.optProb.finalizeDesignVariables()
        self.optProb.finalizeConstraints()
        # Setup initial cache values
        self._setInitialCacheValues()
        self._setSens(sens, sensStep, sensMode, comm)
                
        # We make a split here: If the rank is zero we setup the
        # problem and run SNOPT, otherwise we go to the waiting loop:

        if rank == 0:
            blx, bux, xs = self._assembleContinuousVariables()
            ncon, blc, buc = self._assembleConstraints()
            ff = self._assembleObjective()

            # Initialize the Print and Summary files
            # --------------------------------------
            iPrint = self.options['iPrint'][1]
            PrintFile = os.path.join(self.getOption('outputDirectory'),
                                     self.getOption('Print file'))
            if iPrint != 0:
                ierror = snopt.openunit(iPrint, PrintFile, "replace", "sequential")
                if ierror != 0:
                    raise Error('Failed to properly open %s, ierror = %3d'%
                                (PrintFile,ierror))

            iSumm = self.options['iSumm'][1]
            SummFile = os.path.join(self.getOption('outputDirectory'),
                                    self.getOption('Summary file'))
            if iSumm != 0:
                ierror = snopt.openunit(iSumm, SummFile, "replace", "sequential")
                if ierror != 0:
                    raise Error('Failed to properly open %s, ierror = %3d'%
                                  (SummFile, ierror))
      

            # Determine the sparsity structure of the full jacobian
            # -----------------------------------------------------
            # Get nonlinear part:
            gcon = {}
            for iCon in self.optProb.constraints:
                con = self.optProb.constraints[iCon]
                if not con.linear:
                    gcon[iCon] = con.jac

            fullJacobian = self.optProb.processConstraintJacobian(
                gcon, linearFlag=False)

            # If we have linear constraints those are already done actually.
            if self.optProb.linearJacobian is not None:
                fullJacobian = sparse.vstack([fullJacobian,
                                              self.optProb.linearJacobian])

            # Convert to the csc format that snopt needs
            fullJacobian = fullJacobian.tocsc()

            Acol = fullJacobian.data
            indA = fullJacobian.indices + 1
            locA = fullJacobian.indptr + 1
            
            # Set the number of nonlinear constraints snopt *thinks* we have:
            if self.unconstrained:
                nnCon = 1
            else:
                nnCon = self.optProb.nnCon
            
            # Calculate the length of the work arrays
            # --------------------------------------
            nvar = self.optProb.ndvs
            lencw = 500
            leniw = 500 + 100*(ncon+nvar)
            lenrw = 500 + 200*(ncon+nvar)

            self.options['Total integer workspace'][1] = leniw
            self.options['Total real workspace'][1] = lenrw

            cw = numpy.empty((lencw,8), 'c')
            iw = numpy.zeros(leniw, numpy.intc)
            rw = numpy.zeros(lenrw, numpy.float)
            snopt.sninit(iPrint, iSumm, cw, iw, rw)

            # Memory allocation
            nnObj = nvar
            nnJac = nvar
            iObj = numpy.array(0, numpy.intc)
            neA = len(indA)
            neGcon = neA  # The nonlinear Jacobian and A are the same 
            iExit = 0

            # Set the options into the SNOPT instance
            self._set_snopt_options(iPrint, iSumm, cw, iw, rw)

            mincw, miniw, minrw,cw = snopt.snmemb(iExit, ncon, nvar, 
                                                  neA, neGcon, 
                                                  nnCon, nnJac, nnObj,
                                                  cw, iw, rw)

            if (minrw > lenrw) or (miniw > leniw) or (mincw > lencw):
                if mincw > lencw:
                    lencw = mincw
                    cw = numpy.array((lencw, 8), 'c')
                    tcw[:] = ' '
                if miniw > leniw:
                    leniw = miniw
                    iw = numpy.zeros(leniw, numpy.intc)
                if (minrw > lenrw):
                    lenrw = minrw        
                    rw = numpy.zeros(lenrw, numpy.float)

                snopt.sninit(iPrint, iSumm, cw, iw, rw)

                # snInit resets all the options to the defaults. 
                # Set them again!
                self._set_snopt_options(iPrint, iSumm, cw, iw, rw)  

            # end if

            # Setup argument list values
            start = numpy.array(self.options['Start'][1])
            nName = numpy.array([1], numpy.intc)
            ObjAdd = numpy.array([0.], numpy.float)
            ProbNm = numpy.array(self.optProb.name)        
            xs = numpy.concatenate((xs, numpy.zeros(ncon,numpy.float)))
            bl = numpy.concatenate((blx, blc))
            bu = numpy.concatenate((bux, buc))
            lencu = 1
            leniu = 2
            lenru = 3
            cu = numpy.array(["        "],'c')
            iu = numpy.zeros(leniu, numpy.intc)
            ru = numpy.zeros(lenru, numpy.float)
            hs = numpy.zeros(nvar+ncon, numpy.intc)

            Names = numpy.array(["        "],'c')
            pi = numpy.zeros(ncon, numpy.float)
            rc = numpy.zeros(nvar+ncon, numpy.float)
            inform = numpy.array([-1], numpy.intc)
            mincw = numpy.array([0], numpy.intc)
            miniw = numpy.array([0], numpy.intc)
            minrw = numpy.array([0], numpy.intc)
            nS = numpy.array([0], numpy.intc)
            ninf = numpy.array([0], numpy.intc)
            sinf = numpy.array([0.], numpy.float)

            # Set history
            self._setHistory(storeHistory)

            # Check for warm/cold start --- this is snopt specific so
            # we have a special function for this:
            res1, res2 = self._warmStart(warmStart, coldStart)
            if res1 is not None:
                xs[0:nvar] = res1
            if res2 is not None:
                hs = res2.copy()

            # Setup hot start if necessary
            self._hotStart(storeHistory, hotStart)

            # The snopt c interface
            timeA = time.time()
            snopt.snoptc(start, nnCon, nnObj, nnJac, iObj, ObjAdd, ProbNm,
                         self._userfg_wrap, Acol, indA, locA, bl, bu, 
                         Names, hs, xs, pi, rc, inform, mincw, miniw, minrw, 
                         nS, ninf, sinf, ff, cu, iu, ru, cw, iw, rw)
            optTime = time.time()-timeA
          
            if self.storeHistory:
                # Record the full state of variables, xs and hs such
                # that we could perform a warm start. 
                self.hist.writeData('xs', xs)
                self.hist.writeData('hs', hs)
                self.hist.close()

            if iPrint != 0:
                snopt.closeunit(self.options['iPrint'][1])
            if iSumm != 0:
                snopt.closeunit(self.options['iSumm'][1])

            # Store Results
            sol_inform = {}
            sol_inform['value'] = inform
            sol_inform['text'] = self.informs[inform[0]]

            # Create the optimization solution
            funcEval = 1
            sol = self._createSolution(optTime, funcEval, sol_inform, ff)

            if MPI:
                # Broadcast a -1 to indcate SNOPT has finished
                MPI.COMM_WORLD.bcast(-1, root=0)
            
        else: # We are not on the root process so go into waiting loop:
            self.waitLoop()
            sol = None
        # end if

        # Communicate the solution -- We are back to the point where
        # all processors are back together, so a standard bcast is
        # fine.
        if MPI:
            sol = MPI.COMM_WORLD.bcast(sol)
        
        return  sol

    def _userfg_wrap(self, mode, nnJac, x, fobj, gobj, fcon, gcon, 
                     nState, cu, iu, ru):
        """
        The snopt user function. This is what is actually called from snopt.
        
        Essentially nothing is done in this function, but this funcion
        has to precisely match the signature from fortran so must look
        EXACTLY like this. 

        All we do here is call the generic masterFunc in the baseclass
        which will take care of everything else.
        """

        fail = False
        if mode == 0 or mode == 2:
            fobj, fcon, fail = self.masterFunc(x, ['fobj', 'fcon'])
        if mode == 1:
            if self.getOption('Derivative level') != 0:
                gobj, gcon, fail = self.masterFunc(x, ['gobj', 'gcon'])
        if mode == 2:
            if self.getOption('Derivative level') != 0:
                gobj, gcon, fail2 = self.masterFunc(x, ['gobj', 'gcon'])
                fail = fail or fail2
                
        # Flush the files to the buffer for all the people who like to
        # monitor the residual
        snopt.pyflush(self.getOption('iPrint'))
        snopt.pyflush(self.getOption('iSumm'))

        if fail:
            mode = -1

        return mode, fobj, gobj, fcon, gcon
       
    def _warmStart(self, warmStart, coldStart):
        """
        Internal snopt function to do a warm start or cold start. The
        cold start code resides in the generic base class"""

        xs = None
        hs = None
        if warmStart is not None:
            if os.path.exists(warmStart):
                hist = History(warmStart, flag='r')
                xs_tmp = hist.readData('xs')
                hs_tmp = hist.readData('hs')
                hist.close()
                if xs_tmp is not None and hs_tmp is not None:
                    if len(xs_tmp) == len(xs) and len(hs_tmp) == len(hs):
                        xs = xs_tmp.copy()
                        hs = hs_tmp.copy()
                        # Tell snopt to use this warm start information
                        self.setOption('Start', 'Warm start')
                    else:
                        print('The number of variables or constraints in warmStart file do not \
match the number in the current optimization. Ignorning warmStart file and trying cold start.')
                        coldStart = warmStart

                else:
                    print('No warm start information in file. \'xs\' and \'hs\' must be\
 present in history file. Trying cold start instead.')
                    coldStart = warmStart
            else:
                print('warmStart file not found. Continuing without warm restart')
        # end if (warm start)

        # Now if we have a cold start, we can call the common code:
        if coldStart is not None:
            xs = self._coldStart(coldStart)

        # Return 
        return xs, hs

    def _set_snopt_options(self, iPrint, iSumm, cw, iw, rw):
        """
        Set all the options into SNOPT that have been assigned
        by the user
        """
                
        # Set Options from the local options dictionary
        # ---------------------------------------------
        inform = numpy.array([-1], numpy.intc)
        for item in self.set_options:
            name = item[0]
            value = item[1]

            if isinstance(value, str):
                if (name == 'Start'):
                    if (value == 'Cold'):
                        snopt.snset('Cold start', iPrint, iSumm, inform, 
                                    cw, iw, rw)
                    elif (value == 'Warm'):
                        snopt.snset('Warm start', iPrint, iSumm, inform, 
                                    cw, iw, rw)
                elif (name == 'Problem Type'):
                    if (value == 'Minimize'):
                        snopt.snset('Minimize', iPrint, iSumm, inform, 
                                    cw, iw, rw)
                    elif (value == 'Maximize'):
                        snopt.snset('Maximize', iPrint, iSumm, inform, 
                                    cw, iw, rw)
                    elif (value == 'Feasible point'):
                        snopt.snset('Feasible point', iPrint, iSumm, inform, 
                                    cw, iw, rw)
                elif (name == 'Print file'):
                    snopt.snset(name+' '+'%d'%(iPrint), iPrint, iSumm, inform, 
                                cw, iw, rw)
                elif (name == 'Summary file'):
                    snopt.snset(name+' '+'%d'%(iSumm), iPrint, iSumm, inform, 
                                cw, iw, rw)
                else:
                    snopt.snset(name+' '+value, iPrint, iSumm, inform, 
                                cw, iw, rw)
            elif isinstance(value, float):
                snopt.snsetr(name, value, iPrint, iSumm, inform, cw, iw, rw)
            elif isinstance(value, int):
                snopt.snseti(name, value, iPrint, iSumm, inform, cw, iw, rw)
            elif isinstance(value, type(None)):
                snopt.snset(name, iPrint, iSumm, inform, cw, iw, rw)

        return

    def _on_setOption(self, name, value):
        """
        Set Optimizer Option Value (Optimizer Specific Routine)
        
        Documentation last updated:  May. 07, 2008 - Ruben E. Perez
        """
        
        self.set_options.append([name,value])
        
        
    def _on_getOption(self, name):
        """
        Get Optimizer Option Value (Optimizer Specific Routine)
        
        Documentation last updated:  May. 07, 2008 - Ruben E. Perez
        """
        
        pass
        
    def _on_getInform(self, infocode):
        """
        Get Optimizer Result Information (Optimizer Specific Routine)
        
        Keyword arguments:
        -----------------
        id -> STRING: Option Name
        
        Documentation last updated:  May. 07, 2008 - Ruben E. Perez
        """
        
        # 
        mjr_code = (infocode[0]/10)*10
        mnr_code = infocode[0] - 10*mjr_code
        try:
            inform_text = self.informs[mjr_code]
        except:
            inform_text = 'Unknown Exit Status'
        # end try
        
        return inform_text
        
    def _on_flushFiles(self):
        """
        Flush the Output Files (Optimizer Specific Routine)
        
        Documentation last updated:  August. 09, 2009 - Ruben E. Perez
        """
        
        # 
        iPrint = self.options['iPrint'][1]
        iSumm = self.options['iSumm'][1]
        if (iPrint != 0):
            snopt.pyflush(iPrint)

        if (iSumm != 0):
            snopt.pyflush(iSumm)
 
#==============================================================================
# SNOPT Optimizer Test
#==============================================================================
if __name__ == '__main__':
    
    # Test SNOPT
    print('Testing ...')
    snopt = SNOPT()
    print(snopt)

