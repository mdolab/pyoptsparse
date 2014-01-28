from __future__ import print_function
from __future__ import absolute_import
#!/usr/bin/env python
"""
pyOpt_optimizer

Holds the Python Design Optimization Classes (base and inherited).

Copyright (c) 2008-2013 by pyOpt Developers
All rights reserved.
Revision: 1.1   $Date: 08/05/2008 21:00$


Developers:
-----------
- Dr. Ruben E. Perez (RP)
- Mr. Peter W. Jansen (PJ)

History
-------
    v. 1.0  - Initial Class Creation (RP, 2008)
    v. 1.1  - Pretty Print of Optimization Problems (PJ, 2008)
            - Bug Fix in setOption (PJ, 2008)
"""

__version__ = '$Revision: $'

"""
To Do:
    - deal with g>=0 & g<=0 depending on optimizer !
    - add timing routine in optimizer class ?? better leave @ specific opt level ?
    - parameters class ?
    - How to handle external gradients calculations ?
    - Implement Restart Capability ?
    - self tests ?
    - add Optimizer Info method
"""

from .pyOpt_error import Error

class Optimizer(object):
    
    """
    Abstract Class for Optimizer Object
    """
    
    def __init__(self, name=None, category=None, def_options=None, informs=None, **kwargs):
        
        """
        Optimizer Class Initialization
        
        **Keyword arguments:**
        
        - name -> STR: Optimizer name, *Default* = {}
        - category -> STR: Optimizer category, *Default* = {}
        - def_options -> DICT: Deafult options, *Default* = {}
        - informs -> DICT: Calling routine informations texts, *Default* = {}
        
        Documentation last updated:  Feb. 03, 2011 - Peter W. Jansen
        """
        
        # 
        self.name = name
        self.category = category
        self.options = {}
        self.options['defaults'] = def_options
        self.informs = informs
        
        # Initialize Options
        for key in def_options:
            self.options[key] = def_options[key]

        koptions = kwargs.pop('options', {})
        for key in koptions:
            self.setOption(key, koptions[key])
        
        
    def _on_setOption(self, name, value):
        
        """
        Set Optimizer Option Value (Optimizer Specific Routine)
        
        **Arguments:**
        
        - name -> STR: Option name
        - value ->   : Option value
        
        Documentation last updated:  Feb. 07, 2011 - Peter W. Jansen
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
            raise Error('Received an unknown option: %s'%repr(name))
        
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
    
