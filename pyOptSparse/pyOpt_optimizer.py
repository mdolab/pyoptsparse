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
from .pyOpt_error import Error

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
    
