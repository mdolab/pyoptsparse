#!/usr/bin/env python
'''
pyOpt_gradient

Holds the Python Design Optimization Gradient Calculation Class.

Copyright (c) 2008-2013 by pyOpt Developers
All rights reserved.
Revision: 1.0   $Date: 20/06/2010 21:00$


Developers:
-----------
- Dr. Ruben E. Perez (RP)
- Mr. Peter W. Jansen (PJ)

History
-------
    v. 1.0  - Initial Class Creation (PJ,RP 2010)
'''

__version__ = '$Revision: $'

'''
To Do:
    - add calc fail flag
'''

# =============================================================================
# Standard Python modules
# =============================================================================
import os, sys
import copy
import pdb

# =============================================================================
# External Python modules
# =============================================================================
import numpy

# =============================================================================
# Extension modules
# =============================================================================
#import extension

# =============================================================================
# Misc Definitions
# =============================================================================
eps = 1.0	# define a value for machine precision
while ((eps/2.0 + 1.0) > 1.0):
    eps = eps/2.0
#end
eps = 2.0*eps
#eps = math.ldexp(1,-52)


# =============================================================================
# Gradient Class
# =============================================================================
class Gradient(object):
    
    '''
    Abstract Class for Optimizer Gradient Calculation Object
    '''
    
    def __init__(self, opt_problem, sens_type, sens_mode='', sens_step={}, *args, **kwargs):
        
        '''
        Optimizer Gradient Calculation Class Initialization
        
        **Arguments:**
        
        - opt_problem -> INST: Optimization instance
        - sens_type -> STR/FUNC: Sensitivity type ('FD', 'CS', or function) 
        
        **Keyword arguments:**
        
        - sens_mode -> STR: Parallel flag [''-serial,'pgc'-parallel], *Default* = '' 
        - sens_step -> INT: Step size, *Default* = {} [=1e-6(FD), 1e-20(CS)] 
        
        Documentation last updated:  Feb. 03, 2011 - Peter W. Jansen
        '''
        
        # 
        self.opt_problem = opt_problem
        if isinstance(sens_type,str):
            self.sens_type = sens_type.lower()
        else:
            self.sens_type = sens_type
        #end
        if (sens_step == {}):
            if (self.sens_type == 'fd'):
                self.sens_step = 1.0e-6
            elif (self.sens_type == 'cs'):
                self.sens_step = 1.0e-20
            else:
                self.sens_step = sens_step
            #end
        else:
            self.sens_step = sens_step
        #end
        self.sens_mode = sens_mode.lower()
        
        # MPI Setup
        if (self.sens_mode.lower() == 'pgc'):
            try:
                import mpi4py
                from mpi4py import MPI
            except ImportError:
                print 'Error: mpi4py library failed to import'
            #end
            comm = MPI.COMM_WORLD
            self.nproc = comm.Get_size()
            self.myrank = comm.Get_rank()
            if (mpi4py.__version__[0] == '0'):
                self.Send = comm.Send
                self.Recv = comm.Recv
                self.Bcast = comm.Bcast
            elif (mpi4py.__version__[0] == '1'):
                self.Send = comm.send
                self.Recv = comm.recv
                self.Bcast = comm.bcast
            #end
            self.mydvs = xrange(self.myrank,len(opt_problem._variables.keys()),self.nproc)
        else:
            self.myrank = 0
            self.mydvs = xrange(len(opt_problem._variables.keys()))
        #end
        
        
    def getGrad(self, x, group_ids, f, g, *args, **kwargs):
        
        '''
        Get Gradient
        
        **Arguments:**
        
        - x -> ARRAY: Design variables
        - group_ids -> DICT: Group identifications
        - f -> ARRAY: Objective values
        - g -> ARRAY: Constraint values
        
        Documentation last updated:  Feb. 07, 2011 - Peter W. Jansen
        '''
        
        # 
        opt_problem = self.opt_problem
        sens_type = self.sens_type
        sens_mode = self.sens_mode
        sens_step = self.sens_step
        mydvs = self.mydvs
        myrank = self.myrank
        
        
        # 
        dfi = numpy.zeros([len(opt_problem._objectives.keys()),len(mydvs)],'d')
        dgi = numpy.zeros([len(opt_problem._constraints.keys()),len(mydvs)],'d')
        
        if (sens_type == 'fd'):
            
            # Finite Differences
            dh = sens_step            
            xs = x
            k = 0
            for i in mydvs:
                xh = copy.copy(xs)
                xh[i] += dh
                
                # Variables Groups Handling
                if opt_problem.use_groups:
                    xhg = {}
                    for group in group_ids.keys():
                        if (group_ids[group][1]-group_ids[group][0] == 1):
                            xhg[group] = xh[group_ids[group][0]]
                        else:
                            xhg[group] = xh[group_ids[group][0]:group_ids[group][1]]
                        #end
                    #end
                    xh = xhg
                #end
                
                [fph,gph,fail] = opt_problem.obj_fun(xh, *args, **kwargs)
                if isinstance(fph,float):
                    fph = [fph]
                #end
                
                for j in xrange(len(opt_problem._objectives.keys())):
                    dfi[j,k] = (fph[j] - f[j])/dh
                #end
                for j in xrange(len(opt_problem._constraints.keys())):
                    dgi[j,k] = (gph[j] - g[j])/dh
                #end
                k += 1
            #end
            
        elif (sens_type == 'cs'):
            
            # Complex Step
            cdh = sens_step
            cxs = copy.copy(x)
            k = 0
            for i in mydvs:
                cxh = cxs + numpy.zeros(len(cxs),complex)
                cxh[i] = complex(cxh[i],cdh)
                
                # Variables Groups Handling
                if opt_problem.use_groups:
                    cxhg = {}
                    for group in group_ids.keys():
                        if (group_ids[group][1]-group_ids[group][0] == 1):
                            cxhg[group] = cxh[group_ids[group][0]]
                        else:
                            cxhg[group] = cxh[group_ids[group][0]:group_ids[group][1]]
                        #end
                    #end
                    cxh = cxhg
                #end
                
                [cfph,cgph,fail] = opt_problem.obj_fun(cxh, *args, **kwargs)
                if isinstance(cfph,complex):
                    cfph = [cfph]
                #end
                
                for j in xrange(len(opt_problem._objectives.keys())):
                    dfi[j,k] = cfph[j].imag/cdh
                #end
                for j in xrange(len(opt_problem._constraints.keys())):
                    dgi[j,k] = cgph[j].imag/cdh
                #end
                k += 1
            #end
            
            dfi = dfi.astype(float)
            dgi = dgi.astype(float)
            
        else:
            
            # Variables Groups Handling
            if opt_problem.use_groups:
                xg = {}
                for group in group_ids.keys():
                    if (group_ids[group][1]-group_ids[group][0] == 1):
                        xg[group] = x[group_ids[group][0]]
                    else:
                        xg[group] = x[group_ids[group][0]:group_ids[group][1]]
                    #end
                #end
                xn = xg
            else:
                xn = x
            #end
            
            # User Provided Sensitivities
            [df_user,dg_user,fail] = sens_type(xn, f, g, *args, **kwargs)
            
            if isinstance(df_user,list):
                if len(opt_problem._objectives.keys()) == 1:
                    df_user = [df_user]
                #end
                df_user = numpy.array(df_user)
            #end
            if isinstance(dg_user,list):
                dg_user = numpy.array(dg_user)
            #end
            
            # 
            for i in xrange(len(opt_problem._variables.keys())):
                for j in xrange(len(opt_problem._objectives.keys())):
                    dfi[j,i] = df_user[j,i]
                #end
                for j in xrange(len(opt_problem._constraints.keys())):
                    dgi[j,i] = dg_user[j,i]
                #end
            #end
            
        #end
        
        # MPI Gradient Assembly
        df = numpy.zeros([len(opt_problem._objectives.keys()),len(opt_problem._variables.keys())],'d')
        dg = numpy.zeros([len(opt_problem._constraints.keys()),len(opt_problem._variables.keys())],'d')
        if (sens_mode == 'pgc'):
            if (sens_type == 'fd') or (sens_type == 'cs'):
                if myrank != 0:
                    self.Send([myrank, dfi, dgi],dest=0)
                else:
                    p_results = [[myrank, dfi, dgi]]
                    for proc in xrange(1,self.nproc):
                        p_results.append(self.Recv(source=proc))
                    #end
                #end
                if myrank == 0:
                    for proc in xrange(self.nproc):
                        k = 0
                        for i in xrange(p_results[proc][0],len(opt_problem._variables.keys()),self.nproc):
                            df[:,i] = p_results[proc][1][:,k]
                            dg[:,i] = p_results[proc][2][:,k]
                            k += 1
                        #end
                    #end
                #end
            [df,dg] = self.Bcast([df,dg],root=0)
            #end
        else:
            df = dfi
            dg = dgi
        #end
        
        
        return df,dg
        
        
    def getHess(self, *args, **kwargs):
        
        '''
        Get Hessian
        
        Documentation last updated:  June. 20, 2010 - Ruben E. Perez
        '''
        
        # 
        
        
        # 
        return 
    


#==============================================================================
# Optimizer Gradient Calculation Test
#==============================================================================
if __name__ == '__main__':
    
    # Test Optimizer Gradient Calculation
    print 'Testing Optimizer Gradient Calculation...'
    grd = Gradient()
    
