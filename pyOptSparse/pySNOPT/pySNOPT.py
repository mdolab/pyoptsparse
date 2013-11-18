#! /usr/bin/env python
'''
pySNOPT - A variation of the pySNOPT wrapper specificially designed to
work with sparse optimization problems.

Copyright (c) 2013-2013 by Dr. Gaetan Kenway
All rights reserved.
Revision: 1.0   $Date: 20/09/2013 21:00$


Tested on:
---------
Linux with intel

Developers:
-----------
- Dr. Gaetan Kenway (GKK)
- Dr. Graeme Kennedy (GJK)
History
-------
    v. 0.1    - Initial Wrapper Creation (JM, 2000)
'''

__version__ = '$Revision: $'


# =============================================================================
# SNOPT Library
# =============================================================================
try:
    import snopt
except:
    raise ImportError('SNOPT shared library failed to import')
# end if

# =============================================================================
# Standard Python modules
# =============================================================================
import os
import sys
import copy
import time

# =============================================================================
# External Python modules
# =============================================================================
import numpy
import shelve

# =============================================================================
# Extension modules
# =============================================================================
from pyOptSparse import Optimizer
from pyOptSparse import History
# =============================================================================
# Misc Definitions
# =============================================================================
inf = 10.E+20  # define a value for infinity

# Try to import mpi4py and determine rank
try: 
    from mpi4py import MPI
    rank = MPI.COMM_WORLD.rank
except:
    rank = 0
    MPI = None
# end try

eps = numpy.finfo(1.0).eps
# =============================================================================
# SNOPT Optimizer Class
# =============================================================================
class SNOPT(Optimizer):
    
    '''
    SNOPT Optimizer Class - Inherited from Optimizer Abstract Class
    '''
    
    def __init__(self, *args, **kwargs):
        
        '''
        SNOPT Optimizer Class Initialization
        
        **Keyword arguments:**
        
        Documentation last updated:  Feb. 16, 2010 - Peter W. Jansen
        '''
        
        #
        name = 'SNOPT'
        category = 'Local Optimizer'
        def_opts = {
        # SNOPT Printing Options
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
        Optimizer.__init__(self, name, category, def_opts, informs, *args, **kwargs)

        # The state of the variables and the slacks
        self.x_previous = None
        self.callCounter = 0
        self.start_time = None
        self.time_limit = None

    def __solve__(self, opt_prob, gobj_con=None, store_hst=None, 
                  hot_start=None, warm_start=None, cold_start=None, 
                  time_limit=None, *args, **kwargs):
        
        '''
        Run Optimizer (Optimize Routine)
        - opt_problem -> INST: Optimization instance
        - gobj_con -> FUNC: Gradient function 
        
        **Keyword arguments:**
        
        - store_sol -> BOOL: Store solution in Optimization class flag, *Default* = True 
        - disp_opts -> BOOL: Flag to display options in solution text, *Default* = False
        - store_hst -> STR: Filename to store optimization history, *Default* = None (don't store)
        - hot_start -> STR: Filename to read optimization history and hot start a previous
                            optimization. *Default* = None, don't hot start
        - cold_start -> STR: Filename to read optimization history and do a cold start from a prevous
                            optimization. *Default* = None, don't cold start
        - warm_start -> File to read optimization history and do a warm start from a previous
                            optimization. *Default* = None, don't warm start
        - time_limit -> Number: The time limit in seconds for optimziation. A clean 
                                terimination will be executed ater 'time_limit' seconds.
        Documentation last updated:  Feb. 2, 2011 - Peter W. Jansen
        '''

        # Check that gobj_con is not None if Derivative level is non zero:
        if self.getOption('Derivative level') <> 0 and gobj_con is None:
            print 'Erorr: Derivative level is not 0 and gradient function not supplied'
            sys.exit(0)
        
        # Pull off starting time, if necessary
        self.start_time = time.time()
        if time_limit is not None:
            self.time_limit = time_limit

        # Save the optimization problem and the gradient function
        self.opt_prob = opt_prob
        self.gobj_con = gobj_con

        # We make a split here: If the rank is zero we setup the
        # problem and run SNOPT, otherwise we go to the waiting loop:

        if rank == 0:

            # Get the variable names and variable bounds
            # ------------------------------------------
            blx = []
            bux = []
            xs = []
            for dvSet in self.opt_prob.variables.keys():
                for dvGroup in self.opt_prob.variables[dvSet]:
                    for var in self.opt_prob.variables[dvSet][dvGroup]:
                        if var.type == 'c':
                            blx.append(var.lower)
                            bux.append(var.upper)
                            xs.append(var.value)

                        elif (self.opt_prob._variables[key].type == 'i'):
                            raise IOError('SNOPT cannot handle integer design variables')
                        elif (self.opt_prob._variables[key].type == 'd'):
                            raise IOError('SNOPT cannot handle discrete design variables')
                        # end if
                    # end for
                # end for
            # end for
            blx = numpy.array(blx)
            bux = numpy.array(bux)
            xs = numpy.array(xs)

            # Constraints Handling -- make sure nonlinear constraints go first!
            blc = []
            buc = []
            if len(self.opt_prob.constraints) > 0: 
                for key in self.opt_prob.constraints.keys():
                    if not self.opt_prob.constraints[key].linear:
                        if (self.opt_prob.constraints[key].type == 'i'):
                            blc.extend(self.opt_prob.constraints[key].lower)
                            buc.extend(self.opt_prob.constraints[key].upper)
                        else:
                            print 'Error: only inequality constraints allowed; \
use the same upper/lower bounds for equality constraints'
                            sys.exit(1)
                        # end if
                    # end if
                # end for

                for key in self.opt_prob.constraints.keys():
                    if self.opt_prob.constraints[key].linear:
                        if (self.opt_prob.constraints[key].type == 'i'):
                            blc.extend(self.opt_prob.constraints[key].lower)
                            buc.extend(self.opt_prob.constraints[key].upper)
                        else:
                            print 'Error: only inequality constraints allowed; \
use the same upper/lower bounds for equality constraints'
                            sys.exit(1)
                        # end if
                    # end if
                # end for
            else:
                blc.append(-inf)
                buc.append( inf)
            # end if
            ncon = len(blc)
            blc = numpy.array(blc)
            buc = numpy.array(buc)

           # Objective Handling
            objfunc = self.opt_prob.obj_fun
            nobj = len(self.opt_prob.objectives.keys())
            ff = []
            for key in self.opt_prob.objectives.keys():
                ff.append(self.opt_prob.objectives[key].value)
            #end
            ff = numpy.array(ff)

            # Initialize the Print and Summary files
            # --------------------------------------
            iPrint = self.options['iPrint'][1]
            PrintFile = self.options['Print file'][1]
            if iPrint != 0:
                if os.path.isfile(PrintFile):
                    os.remove(PrintFile)

                ierror = snopt.openunit(iPrint, numpy.array(PrintFile), 
                                        numpy.array('new'), 
                                        numpy.array('sequential'))
                if ierror != 0:
                    raise IOError('Failed to properly open %s, ierror = %3d'%
                                  (PrintFile,ierror))

            iSumm = self.options['iSumm'][1]
            SummFile = self.options['Summary file'][1]
            if iSumm != 0:
                if os.path.isfile(SummFile):
                    os.remove(SummFile)
                ierror = snopt.openunit(iSumm, numpy.array(SummFile), 
                                        numpy.array('new'), 
                                        numpy.array('sequential'))
                if ierror != 0:
                    raise IOError('Failed to properly open %s, ierror = %3d'%
                                  (SummFile,ierror))

            self.opt_prob.reorderConstraintJacobian(reorder=['nonLinear','linear'])

            # We will also assemble just the nonlinear jacobain to
            # determine the number nonzero entries. This must remain
            # fixed in subsequent iterations
            gcon = {}
            for iCon in self.opt_prob.constraints:
                con = self.opt_prob.constraints[iCon]
                if con.dense:
                    gcon[iCon] = con.jac.todense()
                else:
                    gcon[iCon] = con.jac
            # end for
            gobj = numpy.zeros(self.opt_prob.ndvs)

            gobj, fullJacobian = self.opt_prob.processDerivatives(
                gobj, gcon, linearConstraints=True, nonlinearConstraints=True)
            
            fullJacobian = fullJacobian.tocsc()
            Acol = fullJacobian.data
            indA = fullJacobian.indices + 1
            locA = fullJacobian.indptr + 1
            
            self.nnCon = self.opt_prob.nnCon

            # Calculate the length of the work arrays
            # --------------------------------------
            nvar = self.opt_prob.ndvs
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
            nnCon = self.nnCon
            nnJac = nvar
            iObj = numpy.array(0, numpy.intc)
            neA = len(indA)
            neGcon = neA  # The nonlinear Jacobian and A are the same 
            iExit = 0

            # Set the options into the SNOPT instance
            self._set_snopt_options(iPrint, iSumm, cw, iw, rw)

            mincw, miniw, minrw,cw = snopt.snmemb(iExit, ncon, nvar, 
                                                  neA, neGcon, 
                                                  nnCon, nnJac, nnObj, cw, iw, rw)

            if (minrw > lenrw) or (miniw > leniw) or (mincw > lencw):
                print 'pySNOPT: Initial memory estimate for snopt insufficient'
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
            ProbNm = numpy.array(self.opt_prob.name)        
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

            # Open history file if required:
            self.store_hst = False
            if store_hst:
                self.hist = History(store_hst)
                self.store_hst = True
            # end if


            # Check for warm start 
            # ------------------------------------------
            if warm_start is not None:
                if os.path.exists(warm_start):
                    hist = History(warm_start, flag='r')
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
                            print 'The number of variables or constraints in warm_start file do not \
match the number in the current optimization. Ignorning warm_start file and trying cold start.'
                            cold_start = warm_start
                        # end if
                    else:
                        print 'No warm start information in file. \'xs\' and \'hs\' must be\
 present in history file. Trying cold start.'
                        cold_start = warm_start
                else:
                    print 'warm_file not found. Continuing without warm restart'
                # end if
            # end if

            # Check for cold start 
            # ------------------------------------------
            if cold_start is not None:
                if os.path.exists(cold_start):
                    cold_file = shelve.open(cold_start,flag='r')
                    last_key = cold_file['last']
                    x = cold_file[last_key]['x_array'].copy()*self.opt_prob.xscale
                    cold_file.close()
                    if len(x) == nvar:
                        xs[0:nvar] = x.copy()
                    else:
                        print 'The number of variable in cold_start file do not \
match the number in the current optimization. Ignorning cold_start file'
                    # end if
                else:
                    print 'Cold file not found. Continuing without cold restart'
                # end if
            # end if

            self.hot_start = None
            # Determine if we want to do a hot start:
            if hot_start is not None:
                # Now, if if the hot start file and the history are
                # the SAME, we don't allow that. We will create a copy
                # of the hot_start file and use *that* instead. 
                import tempfile, shutil
                if store_hst == hot_start:
                    if os.path.exists(hot_start):
                        fname = tempfile.mktemp()
                        shutil.copyfile(store_hst, fname)
                        self.hot_start = History(fname, temp=True, flag='r')
                else:
                    self.hot_start = History(hot_start, temp=False, flag='r')
                # end if
            # end if

            # The snopt c interface
            snopt.snoptc(start, nnCon, nnObj, nnJac, iObj, ObjAdd, ProbNm,
                         self.userfg_wrap, Acol, indA, locA, bl, bu, 
                         Names, hs, xs, pi, rc, inform, mincw, miniw, minrw, 
                         nS, ninf, sinf, ff, cu, iu, ru, cw, iw, rw)

            if self.store_hst:
                # Record the full state of variables, xs and hs such
                # that we could perform a warm start. 
                self.hist.writeData('xs', xs)
                self.hist.writeData('hs', hs)
                self.hist.close()
            # end if

            if MPI:
                # Broadcast a -1 to indcate SNOPT has finished
                MPI.COMM_WORLD.bcast(-1, root=0)
            # end if

            if iPrint != 0:
                snopt.closeunit(self.options['iPrint'][1])
            if iSumm != 0:
                snopt.closeunit(self.options['iSumm'][1])

            # Store Results
            sol_inform = {}
            sol_inform['value'] = inform
            sol_inform['text'] = self.informs[inform[0]]

        else: # We are not on the rot process so go into waiting loop:

            mode = None
            info = None
            while True:
                # * Note*: No checks for MPI here since this code is
                # * only run in parallel, which assumes mpi4py is working

                # Receive mode and quit if mode is -1:
                mode = MPI.COMM_WORLD.bcast(mode, root=0)
                if mode == -1:
                    break

                # Otherwise receive info from shell function
                info = MPI.COMM_WORLD.bcast(info, root=0)

                # Call the internal function we should have called from
                # SNOPT. We don't care about return values on these procs

                self.userfg(*info)
            # end while

            ff = None
            xs = None
            sol_inform = None
            nvar = None

        # end if

        return 

    def userfg_wrap(self, mode, nnJac, x, fObj, gObj, fCon, gCon, nState, cu, iu, ru):

        '''
        The snopt user function. This is what is actually called from snopt.
        
        It is only called on the root processor actually running
        snopt. History processing is also performed here. The reason
        for this is that the re-runnign on history is serial so it
        makes sense to only read on processor that actually requires
        the data.
        '''
   
        x = x/self.opt_prob.xscale

        # Determine if we've exeeded the time limit:
        if self.time_limit:
            if time.time() - self.start_time > self.time_limit:

                # Broadcast a -1 to indcate SNOPT has finished and
                # make the remainder of the processors finish.
                if MPI:
                    MPI.COMM_WORLD.bcast(-1, root=0)
                # end if
                    
                # Mode = -2 will tell SNOPT that we want to finish
                # (immediately)
                mode = -2
                return mode 
            # end if
        # end if

        # ------------------ Hot Start Processing ------------------
        if self.hot_start:
            if self.hot_start.validPoint(self.callCounter, x):
                data = self.hot_start.read(self.callCounter)
                xn = data['x']
                x_array = data['x_array']
                fObj = data['fobj']
                fCon = data['fcon']
                fail = data['fail']
                if fail: 
                    mode = -1

                fcon_return = self.opt_prob.processConstraints(fCon)

                # Just pass gobj and gcon back if no gradient evaluated
                gobj_return = gObj
                gcon_return = gCon
                gradEvaled = False
                if data.has_key('gobj'):
                    gradEvaled = True
                    gObj = data['gobj']
                    gCon = data['gcon']
                    gobj_return, gcon_return = self.opt_prob.processDerivatives(
                        gObj, gCon, linearConstraints=False, nonlinearConstraints=True)
                    gcon_return = gcon_return.tocsc().data
                # end if

                # Write Data to history (if required):
                if self.store_hst:
                    self.hist.write(self.callCounter, fObj, fCon, fail, xn, x, gradEvaled, 
                                    gObj, gCon, mode=mode, feasibility=ru[0], optimality=ru[1],
                                    merit=ru[1], majorIt=iu[0])

                    
                # end if

                self.callCounter += 1

                return mode, fObj, gobj_return, fcon_return, gcon_return
            # end if

            # We have used up all the information in hot start
            # so we can close the hot start file
            self.hot_start.close()
            self.hot_start = None
        # end if
        # ----------------------------------------------------------

        userfg_args = [mode, nnJac, x, fObj, gObj, fCon, gCon, cu, iu, ru]
        if MPI:
            # Broadcast the type of call (0 means regular call)
            MPI.COMM_WORLD.bcast(0, root=0)

            # Broadcast the requried arguments
            MPI.COMM_WORLD.bcast(userfg_args)
        # end if

        # Call userfg and return result
        return self.userfg(*userfg_args)

    def userfg(self, mode, nnJac, x, fobj, gobj, fcon, gcon, cu, iu, ru):
        '''
        The snopt user function. This is what would normally be called
        from snopt. This function is called on all processors. 
        
        mode = the mode used to indicate whether to eval gradient or not
        njac = the number of Jacobian elements
        x = the design variables
        fobj = objective function
        gobj = gradient of the objective
        fcon = constraints
        gcon = gradient of the constraints
        '''

        xn = self.opt_prob.processX(x)

        # If the gradient isn't calculated, gobj and gcon pass back
        # through
        gobj_return = gobj
        gcon_return = gcon

        # Flush the files to the buffer for all the people who like to 
        # monitor the residual           
        if self.options['iPrint'][1] != 0:
            snopt.pyflush(self.options['iPrint'][1])
        if self.options['iSumm'][1] != 0:
            snopt.pyflush(self.options['iSumm'][1])
            
        # Evaluate the function
        if mode == 0 or mode == 2:
            fargs = self.opt_prob.obj_fun(xn)
            if rank == 0:
                # Process fcon and fobj
                fobj = fargs[0]
                fcon = fargs[1]
                fail = fargs[2]

                fcon_return = self.opt_prob.processConstraints(fcon)
                if fail:
                    mode = -1
                # end if
            # end if
        else:
            if rank == 0:
                fcon_return = fcon.copy()
            # end if
        # end if
            
     
        # Check if last point is *actually* what we evaluated for
        # mode=1
        diff = 1.0
        if not (self.x_previous is None):
            diff = numpy.dot(x-self.x_previous, x-self.x_previous)

        if self.x_previous is None:
            self.x_previous = numpy.zeros(x.size)
        # end if

        self.x_previous[:] = x[:]
        gradEvaled = False
        if self.getOption('Derivative level') <> 0:
            # Evaluate the gradient
            if mode == 2 or (mode == 1 and diff == 0.0):
                # mode == 2: Evaluate the gradient
                # or
                # mode == 1: Only the gradient is required and the previously
                # evaluated point is the same as this point. Evaluate only
                # the gradient                
                gradEvaled = True
                gargs = self.gobj_con(xn, fobj, fcon)

                # Non root rank return for gradient evaluation
                if rank <> 0:
                    return

                # Extract returns
                gobj = gargs[0]
                gcon = gargs[1]
                fail = gargs[2]

                if fail:
                    mode = -1
                # end if

                # Run the standard process derivatives function
                gobj_return, gcon_return = self.opt_prob.processDerivatives(
                    gobj, gcon, linearConstraints=False, nonlinearConstraints=True)
                gcon_return = gcon_return.tocsc().data

            elif mode == 1:
                # mode == 1: only gradient is required, but the
                # previously evaluated point is different. Evaluate the
                # objective then the gradient
                gradEvaled = True
                fargs = self.opt_prob.obj_fun(xn)
                gargs = self.gobj_con(xn, fobj, fcon)

                # Non root rank return for gradient evaluation
                if rank <> 0:
                    return 

                # Extract returns
                gobj = gargs[0]
                gcon = gargs[1]
                fail = gargs[2]

                if fail:
                    mode = -1

                # Run the standard process derivatives function
                gobj_return, gcon_return = self.opt_prob.processDerivatives(
                    gobj, gcon, linearConstraints=False, nonlinearConstraints=True)
                gcon_return = gcon_return.tocsc().data
            # end if
        # end if

        # Non root rank return for objective evaluation only
        if rank <> 0:
            return 

        # Write Data to history:
        if self.store_hst:
            self.hist.write(self.callCounter, fobj, fcon, fail, xn, x, gradEvaled, 
                            gobj, gcon, mode=mode, feasibility=ru[0], optimality=ru[1],
                            merit=ru[1], majorIt=iu[0])
        # end if
        self.callCounter += 1

        return mode, fobj, gobj_return, fcon_return, gcon_return

    def _set_snopt_options(self, iPrint, iSumm, cw, iw, rw):
        '''
        Set all the options into SNOPT that have been assigned
        by the user
        '''
                
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
        '''
        Set Optimizer Option Value (Optimizer Specific Routine)
        
        Documentation last updated:  May. 07, 2008 - Ruben E. Perez
        '''
        
        # 
        self.set_options.append([name,value])
        
        
    def _on_getOption(self, name):
        '''
        Get Optimizer Option Value (Optimizer Specific Routine)
        
        Documentation last updated:  May. 07, 2008 - Ruben E. Perez
        '''
        
        pass
        
    def _on_getInform(self, infocode):
        '''
        Get Optimizer Result Information (Optimizer Specific Routine)
        
        Keyword arguments:
        -----------------
        id -> STRING: Option Name
        
        Documentation last updated:  May. 07, 2008 - Ruben E. Perez
        '''
        
        # 
        mjr_code = (infocode[0]/10)*10
        mnr_code = infocode[0] - 10*mjr_code
        try:
            inform_text = self.informs[mjr_code]
        except:
            inform_text = 'Unknown Exit Status'
        #end
        
        return inform_text
        
    def _on_flushFiles(self):
        '''
        Flush the Output Files (Optimizer Specific Routine)
        
        Documentation last updated:  August. 09, 2009 - Ruben E. Perez
        '''
        
        # 
        iPrint = self.options['iPrint'][1]
        iSumm = self.options['iSumm'][1]
        if (iPrint != 0):
            snopt.pyflush(iPrint)
        #end
        if (iSumm != 0):
            snopt.pyflush(iSumm)
        #end
        
 
#==============================================================================
# SNOPT Optimizer Test
#==============================================================================
if __name__ == '__main__':
    
    # Test SNOPT
    print 'Testing ...'
    snopt = SNOPT()
    print snopt

