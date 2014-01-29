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
import os
import shelve
from .pyOpt_error import Error
from .pyOpt_history import History
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
        
        # Initialize Options
        for key in defOptions:
            self.options[key] = defOptions[key]

        koptions = kwargs.pop('options', {})
        for key in koptions:
            self.setOption(key, koptions[key])

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
            x = coldFile[lastKey]['x_array'].copy()*self.optProb.xscale
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

        # First, apply the master scaling to the design
        # variables. This restores the variables to the range the user
        # is expecting. 
        x = x/self.optProb.xscale

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
                    gon = data['gcon']

                # Process constraints if we have them
                if fcon is not None:
                    fcon = self.optProb.processConstraints(fcon)

                # Process objective gradient
                if gobj is not None:
                    gobj = self.optProb.processObjectiveGradient(gobj)

                # Process constraint gradient
                if gcon is not None:
                    gcon = self.optProb.processConstraintGradient(gcon)

                # Now, gcon is a coo sparse matrix. Depending on what
                # the optimizer wants, we will convert. The
                # conceivable options are: dense (most), csc (snopt), csr
                # (???), or coo (IPOPT)
                if self.jacType == '2ddense':
                    gcon = gcon.todense()
                elif self.jacType == '1ddense':
                    gcon = gcon.todense().flatten()
                elif self.jacType == 'csc':
                    gcon = gcon.tocsc().data
                elif self.jacType == 'csr':
                    gcon = gcon.tocsr().data
                elif self.jacType === 'coo':
                    gcon = gcon.data # Already in coo format
                
                return fobj, fcon, gobj, gcon

                # We can now safely increment the call counter
                self.callCounter += 1



            # end if
        




        pass


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
    
