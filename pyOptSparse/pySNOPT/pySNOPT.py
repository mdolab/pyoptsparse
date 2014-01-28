from __future__ import absolute_import
from __future__ import print_function
#/bin/env python
'''
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
'''
# =============================================================================
# SNOPT Library
# =============================================================================
try:
    import snopt
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
from ..pyOpt_gradient import Gradient
from ..pyOpt_solution import Solution
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
# end try

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
        '''
        
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
        self.startTime = None
        self.timeLimit = None

    def __call__(self, optProb, sens=None, sensStep=None, sensMode=None,
                  storeHistory=None, hotStart=None, warmStart=None,
                  coldStart=None, timeLimit=None, comm=None):
        '''
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
            '''
        
        self.callCounter = 0
        # Pull off starting time, if necessary
        self.startTime = time.time()
        if timeLimit is not None:
            self.timeLimit = timeLimit

        self.unconstrained = False
        if len(optProb.constraints) == 0:
            # If the user *actually* has an unconstrained problem,
            # snopt sort of chokes with that....it has to have at
            # least one constraint. So we will add one
            # automatically here:
            self.unconstrained = True
            
        # Save the optimization problem and finialize constraint
        # jacobian, in general can only do on root proc
        self.optProb = optProb
        self.optProb._finalizeDesignVariables()
        self.optProb.reorderConstraintJacobian(
            reorder=['nonLinear','linear'])

        # Next we determine what to what to do about
        # derivatives. SNOPT is a little special actually since it can
        # do derivatives itself.
        if sens is None:
            # Tell snopt it should do the derivatives
            self.setOption('Derivative level', 0)
            self.sens = None
        elif isinstance(sens, types.FunctionType):
            # We have function handle for gradients! Excellent!
            self.sens = sens
        elif sens.lower() in ['fd','cs']:
            # Create the gradient class that will operate just like if
            # the user supplied fucntion
            self.sens = Gradient(optProb, sens.lower(), sensStep,
                                 sensMode, comm)
        else:
            raise Error('Unknown value given for sens. Must be None, \'FD\', \
            \'CS\' or a python function handle')
        # end if
                
        # We make a split here: If the rank is zero we setup the
        # problem and run SNOPT, otherwise we go to the waiting loop:

        if rank == 0:

            # Get the variable names and variable bounds
            # ------------------------------------------
            blx = []
            bux = []
            xs = []
            for dvSet in self.optProb.variables.keys():
                for dvGroup in self.optProb.variables[dvSet]:
                    for var in self.optProb.variables[dvSet][dvGroup]:
                        if var.type == 'c':
                            blx.append(var.lower)
                            bux.append(var.upper)
                            xs.append(var.value)

                        elif (self.optProb.variables[key].type == 'i'):
                            raise IOError('SNOPT cannot handle integer design variables')
                        elif (self.optProb.variables[key].type == 'd'):
                            raise IOError('SNOPT cannot handle discrete design variables')
                        # end if
                    # end for
                # end for
            # end for
            blx = numpy.array(blx)
            bux = numpy.array(bux)
            xs = numpy.array(xs)
            
            # Constraints Handling -- make sure nonlinear constraints
            # go first -- this is particular to snopt
            blc = []
            buc = []
                
            for key in self.optProb.constraints.keys():
                if not self.optProb.constraints[key].linear:
                    blc.extend(self.optProb.constraints[key].lower)
                    buc.extend(self.optProb.constraints[key].upper)
         
            for key in self.optProb.constraints.keys():
                if self.optProb.constraints[key].linear:
                    blc.extend(self.optProb.constraints[key].lower)
                    buc.extend(self.optProb.constraints[key].upper)

            if self.unconstrained:
                blc.append(-inf)
                buc.append(inf)
            
            ncon = len(blc)
            blc = numpy.array(blc)
            buc = numpy.array(buc)

            # Objective Handling
            objfunc = self.optProb.objFun
            nobj = len(self.optProb.objectives.keys())
            ff = []
            for key in self.optProb.objectives.keys():
                ff.append(self.optProb.objectives[key].value)

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

      
            # We will also assemble just the nonlinear jacobain to
            # determine the number nonzero entries. This must remain
            # fixed in subsequent iterations
            gcon = {}
            for iCon in self.optProb.constraints:
                con = self.optProb.constraints[iCon]
                gcon[iCon] = con.jac
            # end for
            gobj = numpy.zeros(self.optProb.ndvs)

            if not self.unconstrained:
                gobj, fullJacobian = self.optProb.processDerivatives(
                    gobj, gcon, linearConstraints=True, nonlinearConstraints=True)

                fullJacobian = fullJacobian.tocsc()
                self.nnCon = self.optProb.nnCon
            else:
                fullJacobian = sparse.csc_matrix(numpy.ones(self.optProb.ndvs))
                self.nnCon = 1
            # end if
            
            Acol = fullJacobian.data
            indA = fullJacobian.indices + 1
            locA = fullJacobian.indptr + 1
                
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

            # Open history file if required:
            self.storeHistory = False
            if storeHistory:
                self.hist = History(storeHistory)
                self.storeHistory = True
            # end if

            # Check for warm start 
            # ------------------------------------------
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
                        # end if
                    else:
                        print('No warm start information in file. \'xs\' and \'hs\' must be\
 present in history file. Trying cold start.')
                        coldStart = warmStart
                else:
                    print('warm_file not found. Continuing without warm restart')
                # end if
            # end if

            # Check for cold start 
            # ------------------------------------------
            if coldStart is not None:
                if os.path.exists(coldStart):
                    cold_file = shelve.open(coldStart,flag='r')
                    last_key = cold_file['last']
                    x = cold_file[last_key]['x_array'].copy()*self.optProb.xscale
                    cold_file.close()
                    if len(x) == nvar:
                        xs[0:nvar] = x.copy()
                    else:
                        print('The number of variable in coldStart file do not \
match the number in the current optimization. Ignorning coldStart file')
                    # end if
                else:
                    print('Cold file not found. Continuing without cold restart')
                # end if
            # end if

            self.hotStart = None
     
            # Determine if we want to do a hot start:
            if hotStart is not None:
                # Now, if if the hot start file and the history are
                # the SAME, we don't allow that. We will create a copy
                # of the hotStart file and use *that* instead. 
                import tempfile, shutil
                if storeHistory == hotStart:
                    if os.path.exists(hotStart):
                        fname = tempfile.mktemp()
                        shutil.copyfile(storeHistory, fname)
                        self.hotStart = History(fname, temp=True, flag='r')
                else:
                    self.hotStart = History(hotStart, temp=False, flag='r')
                # end if
            # end if

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

            # Create the optimization solution
            sol = Solution(self.optProb, optTime, 1, sol_inform)

            # Now set the x-values:
            i = 0
            for dvSet in sol.variables.keys():
                for dvGroup in sol.variables[dvSet]:
                    for var in sol.variables[dvSet][dvGroup]:
                        var.value = xs[i]
                        i += 1
            sol.fStar = ff
        else: # We are not on the root process so go into waiting loop:

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
            sol = None
        # end if

        # Communicate the solution
        if MPI:
            sol = MPI.COMM_WORLD.bcast(sol)
        
        return  sol

    def _userfg_wrap(self, mode, nnJac, x, fObj, gObj, fCon, gCon, nState, cu, iu, ru):
        '''
        The snopt user function. This is what is actually called from snopt.
        
        It is only called on the root processor actually running
        snopt. History processing is also performed here. The reason
        for this is that the re-runnign on history is serial so it
        makes sense to only read on processor that actually requires
        the data.
        '''
        
        x = x/self.optProb.xscale

        # Determine if we've exeeded the time limit:
        if self.timeLimit:
            if time.time() - self.startTime > self.timeLimit:

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
        if self.hotStart:
            if self.hotStart.validPoint(self.callCounter, x):
                data = self.hotStart.read(self.callCounter)
                xn = data['x']
                x_array = data['x_array']
                fObj = data['fobj']
                fCon = data['fcon']
                fail = data['fail']
                if fail: 
                    mode = -1

                if not self.unconstrained:
                    fcon_return = self.optProb.processConstraints(fCon)
                else:
                    fcon_return = [0]
        
                # Just pass gobj and gcon back if no gradient evaluated
                gobj_return = gObj
                gcon_return = gCon
                gradEvaled = False
                if data.has_key('gobj'):
                    gradEvaled = True
                    gObj = data['gobj']
                    gCon = data['gcon']
                    gobj_return, gcon_return = self.optProb.processDerivatives(
                        gObj, gCon, linearConstraints=False, nonlinearConstraints=True)
                    if self.unconstrained:
                        gcon_return = numpy.zeros(self.optProb.ndvs)
                    else:
                        gcon_return = gcon_return.tocsc().data
                # end if

                # Write Data to history (if required):
                if self.storeHistory:
                    self.hist.write(self.callCounter, fObj, fCon, fail, xn, x, gradEvaled, 
                                    gObj, gCon, mode=mode, feasibility=ru[0], optimality=ru[1],
                                    merit=ru[1], majorIt=iu[0])
                # end if

                self.callCounter += 1

                return mode, fObj, gobj_return, fcon_return, gcon_return
            # end if

            # We have used up all the information in hot start
            # so we can close the hot start file
            self.hotStart.close()
            self.hotStart = None
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

        xn = self.optProb.processX(x)

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

        # Assume fail is False. This is really subtle: When using
        # snopt with its own internal FD calc, even with derivative
        # level 0 and not changing the gobj, and gcon arrays, snopt
        # will **STILL** call with mode 1, just because. When mode 1
        # is called, with derivative level 0, nothing happens, (no
        # function or gradient is called since this is call is realy
        # meaningless --- compute gradient when no gradient is
        # supplied) and the fail flag isn't set anywhere. That's why it
        # needs to be set here.
        fail=False
        
        # Evaluate the function
        if mode == 0 or mode == 2:
            fargs = self.optProb.objFun(xn)
            if rank == 0:
                # Process fcon and fobj
                fobj = fargs[0]
                fcon = fargs[1]
                fail = fargs[2]
                if not self.unconstrained:
                    fcon_return = self.optProb.processConstraints(fcon)
                else:
                    fcon_return = [0]
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
        if self.getOption('Derivative level') != 0:
            # Evaluate the gradient
            if mode == 2 or (mode == 1 and diff == 0.0):
                # mode == 2: Evaluate the gradient
                # or
                # mode == 1: Only the gradient is required and the previously
                # evaluated point is the same as this point. Evaluate only
                # the gradient                
                gradEvaled = True
                gargs = self.sens(xn, fobj, fcon)

                # Non root rank return for gradient evaluation
                if rank != 0:
                    return

                # Extract returns
                gobj = gargs[0]
                gcon = gargs[1]
                fail = gargs[2]

                if fail:
                    mode = -1
                # end if

                # Run the standard process derivatives function
                gobj_return, gcon_return = self.optProb.processDerivatives(
                    gobj, gcon, linearConstraints=False, nonlinearConstraints=True)

                if self.unconstrained:
                    gcon_return = numpy.zeros(self.optProb.ndvs)
                else:
                    gcon_return = gcon_return.tocsc().data

            elif mode == 1:
                # mode == 1: only gradient is required, but the
                # previously evaluated point is different. Evaluate the
                # objective then the gradient
                gradEvaled = True
                fargs = self.optProb.objFun(xn)
                gargs = self.sens(xn, fobj, fcon)

                # Non root rank return for gradient evaluation
                if rank != 0:
                    return 

                # Extract returns
                gobj = gargs[0]
                gcon = gargs[1]
                fail = gargs[2]

                if fail:
                    mode = -1

                # Run the standard process derivatives function
                gobj_return, gcon_return = self.optProb.processDerivatives(
                    gobj, gcon, linearConstraints=False, nonlinearConstraints=True)

                if self.unconstrained:
                    gcon_return = numpy.zeros(self.optProb.ndvs)
                else:
                    gcon_return = gcon_return.tocsc().data
            # end if
        # end if

        # Non root rank return for objective evaluation only
        if rank != 0:
            return 

        # Write Data to history:
        if self.storeHistory:
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
        # end try
        
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

