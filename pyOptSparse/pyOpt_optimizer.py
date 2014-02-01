#!/usr/bin/env python
"""
pyOpt_optimizer

Holds the Python Design Optimization Classes (base and inherited).

Copyright (c) 2008-2013 by pyOpt Developers
All rights reserved.
Revision: 1.1   $Date: 08/05/2008 21:00$

Developers:
-----------
- Dr. Gaetan K.W. Kenway (GKK)
"""
from __future__ import print_function
# =============================================================================
# Imports 
# =============================================================================
import os, sys, copy, types
import shelve
import numpy
from scipy import sparse
from .pyOpt_gradient import Gradient
from .pyOpt_error import Error
from .pyOpt_history import History
eps = numpy.finfo(1.0).eps

# Try to import mpi4py and determine rank
try: 
    from mpi4py import MPI
    rank = MPI.COMM_WORLD.rank
except:
    rank = 0
    MPI = None
# end try

# =============================================================================
# Optimizer Class
# =============================================================================
class Optimizer(object):
    """
    Base optimizer class

    Parameters
    ----------
    name : str
        Optimizer name
    category : str
        Typicaly local or gobal
    defOptions : dictionary
        A dictionary containing the default options
    informs : dict
        Dictionary of the inform codes
        """
    def __init__(self, name=None, category=None, defOptions=None,
                 informs=None, **kwargs):

        self.name = name
        self.category = category
        self.options = {}
        self.options['defaults'] = defOptions
        self.informs = informs
        self.callCounter = 0
        self.sens = None
        # Initialize Options
        for key in defOptions:
            self.options[key] = defOptions[key]

        koptions = kwargs.pop('options', {})
        for key in koptions:
            self.setOption(key, koptions[key])

        self.optProb = None
        # Default options:
        self.appendLinearConstraints = False
        self.jacType = 'dense'

        # Cache storage 
        self.cache = {'x':None, 'fobj':None, 'fcon':None,
                      'gobj':None, 'gcon':None}
        
    def _setSens(self, sens, sensStep, sensMode, comm):
        """
        Common function to setup sens function
        """

        # Next we determine what to what to do about
        # derivatives. We must hvae a function or we use FD or CS:
        if sens is None:
            if self.name in ['SNOPT']:
                # SNOPT is the only one where None is ok. 
                self.setOption('Derivative level', 0)
                self.sens = None
            else:
                raise Error('\'None\' value given for sens. Must be one \
of \'FD\' or \'CS\' or a user supplied function.')
        elif isinstance(sens, types.FunctionType):
            # We have function handle for gradients! Excellent!
            self.sens = sens
        elif sens.lower() in ['fd','cs']:
            # Create the gradient class that will operate just like if
            # the user supplied fucntion
            self.sens = Gradient(self.optProb, sens.lower(), sensStep,
                                 sensMode, comm)
        else:
            raise Error('Unknown value given for sens. Must be None, \'FD\', \
            \'CS\' or a python function handle')
        # end if
 

    def _coldStart(self, coldStart):
        """
        Common code to do cold restarting. 

        Parameters
        ----------
        coldFile : str
           Filename of the history file to use for the cold start
           """
        xCold = None
        if os.path.exists(coldStart):
            # Note we open in read only mode just in case. We don't
            # have to write anyway
            coldFile = shelve.open(coldStart, flag='r')
            lastKey = coldFile['last']
            x = coldFile[lastKey]['x']
            coldFile.close()
            if len(x) == self.optProb.ndvs:
                xCold = x.copy()
            else:
                print('The number of variable in coldStart file do not \
match the number in the current optimization. Ignorning coldStart file')
                    # end if
        else:
            print('Cold restart file not found. Continuing without cold restart')

        return xCold
        
    def _hotStart(self, storeHistory, hotStart):
        """
        Generic routine for setting up the hot start information
        
        Parameters
        ----------
        storeHistory : str
            File for possible history file. Or None if not writing file.

        hotStart : str
            Filename for history file for hot start
            """

        # By default no hot start
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

    def _setHistory(self, storeHistory):
        """Generic routine for setting history file if required
        
        Parameters
        ----------
        storeHistory : str
            Filename for history file for this optimization
            """
        
        self.storeHistory = False
        if storeHistory:
            self.hist = History(storeHistory)
            self.storeHistory = True

    def masterFunc(self, x, evaluate, **kwargs):
        """
        This is the master function that **ALL** optimizers call from
        the specific signature functions. The reason for this is that
        we can generically do the hot-start replay, history storage,
        timing and possibly chaching once for all optimizers. It also
        takes care of the MPI communcation that allows the optimizer
        to run on one process only, but within a larger MPI context. 

        It does add one additional level of call, but we think it is
        well worth it for reduce code duplication
        
        Parameters
        ----------
        x : array
            This is the raw x-array data from the optimizer
        evaluate : list of strings
            This list containts at least one of 'fobj', 'fcon', 'gobj'
            or 'gcon'. This list tells this function which of the
            values is required on return
            """

        # We are hot starting, we should be able to read the required
        # information out of the hot start file, process it and then
        # fire it back to the specific optimizer
        
        if self.hotStart:
            if self.hotStart.validPoint(self.callCounter, x):
                data = self.hotStart.read(self.callCounter)

                if self.storeHistory:
                    # Just dump the (exact) dictionary back out:
                    self.hist.write(self.callCounter, data)

                # Since we know it is a valid point, we can be sure
                # that it contains the information we need
                fobj = None
                fcon = None
                gobj = None
                gcon = None

                if 'fobj' in evaluate:
                    fobj = data['fobj']
                if 'fcon' in evaluate:
                    fcon = data['fcon']
                if 'gobj' in evaluate:
                    gobj = data['gobj']
                if 'gcon' in evaluate:
                    gcon = data['gcon']
                fail = data['fail']
                
                returns = []
                
                # Process objective if we have one (them)
                if fobj is not None:
                    returns.append(self.optProb.processObjective(fobj))
                    
                # Process constraints if we have them
                if fcon is not None:
                    fcon = self.optProb.processNonlinearConstraints(fcon)
                    if self.appendLinearConstraints:
                        fcon = numpy.append(fcon, self.optProb.evaluateLinearConstraints(x))
                    returns.append(fcon)
                    
                # Process objective gradient
                if gobj is not None:
                    returns.append(self.optProb.processObjectiveGradient(gobj))
                    
                # Process constraint gradient
                if gcon is not None:
                    gcon = self.optProb.processConstraintJacobian(gcon)
                    if self.appendLinearConstraints:
                        gcon = sparse.vstack([gcon,
                                              self.optProb.linearJacobian])
                    gcon = self.convertJacobian(gcon)
                    returns.append(gcon)
                    
                # We can now safely increment the call counter
                self.callCounter += 1
                returns.append(fail)

                return returns
            # end if (valid hot start point)  

            # We have used up all the information in hot start so we
            # can close the hot start file
            self.hotStart.close()
            self.hotStart = None
        # end if (hot starting)

        # Now we have to actually run our function...this is where the
        # MPI gets a little tricy. Up until now, only the root proc
        # has called up to here...the rest of them are waiting at a
        # broadcast to know what to do. 

        args = [x, evaluate]
        if MPI:
            # Broadcast the type of call (0 means regular call)
            MPI.COMM_WORLD.bcast(0, root=0)

            # Now broadcast out the required arguments:
            MPI.COMM_WORLD.bcast(args)
        # end if

        return self.masterFunc2(*args)

    def masterFunc2(self, x, evaluate, writeHist=True):
        """
        Another shell function. This function is now actually called
        on all the processors.
        """

        # Our goal in this function is to return the values requested
        # in 'evaluate' for the corresponding x. We have to be a
        # little cheaky here since some optimizers will make multiple
        # call backs with the same x, one for the objective and one
        # for the constraint. We therefore at the end of each function
        # or sensitivity call we cache the x value and the fobj, fcon,
        # gobj, and gcon values such that on the next pass we can just
        # read them and return.

        xscaled = x/self.optProb.xscale
        xuser = self.optProb.processX(xscaled)

        masterFail = False

        # Set basic parameters in history
        hist = {'x':x, 'xuser':xuser, 'xscaled':xscaled}
        returns = []
        # Start with fobj:

        if 'fobj' in evaluate:
            if numpy.linalg.norm(x-self.cache['x']) > eps:
                fobj, fcon, fail = self.optProb.objFun(xuser)
                # User values stored is immediately 
                self.cache['fobj_user'] = copy.deepcopy(fobj)
                self.cache['fcon_user'] = copy.deepcopy(fcon)
                
                # Process constraints 
                fcon = self.optProb.processNonlinearConstraints(fcon)
                if self.appendLinearConstraints:
                    fcon = numpy.append(fcon, self.optProb.evaluateLinearConstraints(x))

                # Now clear out gobj and gcon in the cache since these
                # are out of date and set the current ones
                self.cache['gobj'] = None
                self.cache['gcon'] = None
                self.cache['x'] = x.copy()
                self.cache['fobj'] = copy.deepcopy(fobj)
                self.cache['fcon'] = copy.deepcopy(fcon)

                # Update fail flag
                masterFail = masterFail or fail

            # fobj is now in cache
            returns.append(self.cache['fobj'])
            hist['fobj'] = self.cache['fobj_user']
            
        if 'fcon' in evaluate:
            if numpy.linalg.norm(x-self.cache['x']) > eps:
                fobj, fcon, fail = self.optProb.objFun(xuser)
                # User values stored is immediately 
                self.cache['fobj_user'] = copy.deepcopy(fobj)
                self.cache['fcon_user'] = copy.deepcopy(fcon)

                # Process constraints
                fcon = self.optProb.processNonlinearConstraints(fcon)
                if self.appendLinearConstraints:
                    fcon = numpy.append(fcon, self.optProb.evaluateLinearConstraints(x))

                # Now clear out gobj and gcon in the cache since these
                # are out of date and set the current ones
                self.cache['gobj'] = None
                self.cache['gcon'] = None
                self.cache['x'] = x.copy()
                self.cache['fobj'] = copy.deepcopy(fobj)
                self.cache['fcon'] = copy.deepcopy(fcon)

                # Update fail flag
                masterFail = masterFail or fail

            # fcon is now in cache
            returns.append(self.cache['fcon'])
            hist['fcon'] = self.cache['fcon_user']
            
        if 'gobj' in evaluate:
            if numpy.linalg.norm(x-self.cache['x']) > eps:
                # Previous evaluated point is *different* than the
                # point requested for the derivative. Recusively call
                # the routine with ['fobj', and 'fcon']
                self.masterFunc2(x, ['fobj', 'fcon'], writeHist=False)

            # Now, the point has been evaluated correctly so we
            # determine if we have to run the sens calc:

            if self.cache['gobj'] is None:
                gobj, gcon, fail = self.sens(
                    xuser, self.cache['fobj'], self.cache['fcon_user'])
                # User values are stored is immediately 
                self.cache['gobj_user'] = copy.deepcopy(gobj)
                self.cache['gcon_user'] = copy.deepcopy(gcon)

                # Process objective gradient for optimizer
                gobj = self.optProb.processObjectiveGradient(gobj)
         
                # Process constraint gradients for optimizer
                gcon = self.optProb.processConstraintJacobian(gcon)
                if self.appendLinearConstraints:
                    gcon = sparse.vstack([gcon,
                                          self.optProb.linearJacobian])
                gcon = self.convertJacobian(gcon)
                # Set the cache values:
                self.cache['gobj'] = gobj.copy()
                self.cache['gcon'] = gcon.copy()

                # Update fail flag
                masterFail = masterFail or fail
                
            # gobj is now in the cache
            returns.append(self.cache['gobj'])
            hist['gobj'] = self.cache['gobj_user']
                
        if 'gcon' in evaluate:
            if numpy.linalg.norm(x-self.cache['x']) > eps:
                # Previous evaluated point is *different* than the
                # point requested for the derivative. Recusively call
                # the routine with ['fobj', and 'fcon']
                self.masterFunc2(x, ['fobj', 'fcon'], writeHist=False)
         
            # Now, the point has been evaluated correctly so we
            # determine if we have to run the sens calc:
            if self.cache['gcon'] is None:
                gobj, gcon, fail = self.sens(
                    xuser, self.cache['fobj'], self.cache['fcon_user'])
                # User values stored is immediately 
                self.cache['gobj_user'] = copy.deepcopy(gobj)
                self.cache['gcon_user'] = copy.deepcopy(gcon)

                # Process objective gradient for optimizer
                gobj = self.optProb.processObjectiveGradient(gobj)

                # Process constraint gradients for optimizer
                gcon = self.optProb.processConstraintJacobian(gcon)
                if self.appendLinearConstraints:
                    gcon = sparse.vstack([gcon,
                                          self.optProb.linearJacobian])
                gcon = self.convertJacobian(gcon)

                # Set cache values
                self.cache['gobj'] = gobj.copy()
                self.cache['gcon'] = gcon.copy()

                # Update fail flag
                masterFail = masterFail or fail

            # gcon is now in the cache
            returns.append(self.cache['gcon'])
            hist['gcon'] = self.cache['gcon_user']
                
        # Put the fail flag in the history:
        hist['fail'] = masterFail

        # Write history if necessary
        if rank == 0 and writeHist and self.storeHistory:
            self.hist.write(self.callCounter, hist)

        # We can now safely increment the call counter
        self.callCounter += 1

        # Tack the fail flag on at the end
        returns.append(masterFail)

        return returns

    def convertJacobian(self, gcon):
        """
        Convert gcon which is a coo matrix into the format we need
        """

        # Now, gcon is a coo sparse matrix. Depending on what the
        # optimizer wants, we will convert. The conceivable options
        # are: dense (most), csc (snopt), csr (???), or coo (IPOPT)
        if self.jacType == '2ddense':
            gcon = gcon.todense()
        elif self.jacType == '1ddense':
            gcon = gcon.todense().flatten()
        elif self.jacType == 'csc':
            gcon = gcon.tocsc().data
        elif self.jacType == 'csr':
            gcon = gcon.tocsr().data
        elif self.jacType == 'coo':
            gcon = gcon.data # Already in coo format

        return gcon
            
    def waitLoop(self):
        """Non-root processors go into this waiting loop while the
        root proc does all the work in the optimization algorithm"""
        
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

            # Call the generic internal function. We don't care
            # about return values on these procs
            self.masterFunc2(*info)

    def _setInitialCacheValues(self):
        """
        Once we know that the optProb has been set, we populate the
        cache with a magic numbers. If the starting points for your
        optimization is -9999999999 then you out of luck!
        """
        self.cache['x'] = -999999999*numpy.ones(self.optProb.ndvs)

    def _assembleContinuousVariables(self):
        """
        Utility function for assembling the design variables. Most
        optimizers here use continuous variables so this chunk of code
        can be reused. 
        """

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
                        raise Error('%s cannot handle integer design variables'% self.name)
                    elif (self.optProb.variables[key].type == 'd'):
                        raise Error('%s cannot handle discrete design variables'% self.name)

        blx = numpy.array(blx)
        bux = numpy.array(bux)
        xs = numpy.array(xs)
        
        return blx, bux, xs

    def _assembleConstraints(self):
        """
        Utility function for assembling the design variables. Most
        optimizers here use continuous variables so this chunk of code
        can be reused. 
        """

        # Constraints Handling -- make sure nonlinear constraints
        # go first -- this is particular to slsqp
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

        return ncon, blc, buc

    def _assembleObjective(self):
        """
        Utility function for assembling the design variables. Most
        optimizers here use continuous variables so this chunk of code
        can be reused. 
        """

        nobj = len(self.optProb.objectives.keys())
        if nobj == 0:
            # NO objective, add one                                                                            
            self.optProb.addObj('f', scale=1.0)
        elif nobj <> 1:
            raise Error('%s can only use one objective'% self.name)
        ff = numpy.array([0.0])

        return ff

    def _createSolution(self, optTime, funcEval, sol_inform, obj):
        """
        Generic routine to create the solution after an optimizer
        finishes.
        """
        sol = Solution(self.optProb, optTime, funcEval, sol_inform)

        # Now set the x-values:
        i = 0
        for dvSet in sol.variables.keys():
            for dvGroup in sol.variables[dvSet]:
                for var in sol.variables[dvSet][dvGroup]:
                    var.value = xs[i]
                    i += 1
        sol.fStar = obj

        return sol
        
    def _on_setOption(self, name, value):
        """
        Set Optimizer Option Value (Optimizer Specific Routine)
        """
        raise Error('This optimizer hsa not implemented _on_setOption')
        
    def setOption(self, name, value=None):
        """
        Generic routine for all option setting. This routine does
        error checking on the type of the value. 

        Parameters
        ----------
        name : str
            Name of the option to set
        value : varies
            Variable value to set. 
            """

        if name in self.options['defaults']:
            if type(value) == self.options['defaults'][name][0]:
                self.options[name] = [type(value), value]
            else:
                raise Error('Value type for option %s was incorrect. It was \
expecting type \'%s\' by received type \'%s\''% (
                        name, self.options['defaults'][name][0], type(value)))
        else:
            raise Error('Received an unknown option: %s'% repr(name))
        
        # Now call the optimizer specific routine
        self._on_setOption(name, value)
        
    def _on_getOption(self, name):
        """
        Routine to be implemented by optimizer
        """
        raise Error('This optimizer haa not implemented _on_getOption')
        
    def getOption(self, name):
        """
        Return the optimizer option value for name

        Parameters
        ----------
        name : str
            name of option for which to retrieve value

        Returns
        -------
        value : varies
            value of option for 'name'
            """
        
        if name in self.options['defaults']:
            return self.options[name][1]
        else:	
            raise Error('Received an unknown option: %s.'%repr(name))

        # Now call the optimizer specific routine
        self._on_getOption(name)
        
    def _on_getInform(self, info):
        """
        Routine to be implemented by optimizer
        """        
        raise Error('This optimizer has not implemented _on_getInform')
        
    def getInform(self, infocode=None):
        """
        Get optimizer result infom code at exit

        Parameters
        ----------
        infocode : int
            Integer information code
            """

        if infocode is None:
            return self.informs
        else:
            return self._on_getInform(infocode)
        
#==============================================================================
# Optimizer Test
#==============================================================================
if __name__ == '__main__':
    
    # Test Optimizer
    print('Testing Optimizer...')
    opt = Optimizer()
    
