#!/usr/bin/env python
'''
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
'''

__version__ = '$Revision: $'

'''
To Do:
    - deal with g>=0 & g<=0 depending on optimizer !
    - add timing routine in optimizer class ?? better leave @ specific opt level ?
    - parameters class ?
    - How to handle external gradients calculations ?
    - Implement Restart Capability ?
    - self tests ?
    - add Optimizer Info method
'''

# =============================================================================
# Standard Python modules
# =============================================================================
import os, sys

# =============================================================================
# External Python modules
# =============================================================================
#import external

# =============================================================================
# Extension modules
# =============================================================================
from pyOpt_optimization import Optimization

# =============================================================================
# Optimizer Class
# =============================================================================
class Optimizer(object):
    
    '''
    Abstract Class for Optimizer Object
    '''
    
    def __init__(self, name={}, category={}, def_options={}, informs={}, *args, **kwargs):
        
        '''
        Optimizer Class Initialization
        
        **Keyword arguments:**
        
        - name -> STR: Optimizer name, *Default* = {}
        - category -> STR: Optimizer category, *Default* = {}
        - def_options -> DICT: Deafult options, *Default* = {}
        - informs -> DICT: Calling routine informations texts, *Default* = {}
        
        Documentation last updated:  Feb. 03, 2011 - Peter W. Jansen
        '''
        
        # 
        self.name = name
        self.category = category
        self.options = {}
        self.options['defaults'] = def_options
        self.informs = informs
        
        # Initialize Options
        def_keys = def_options.keys()
        for key in def_keys:
            self.options[key] = def_options[key]
        #end
        koptions = kwargs.pop('options',{})
        kopt_keys = koptions.keys()
        for key in kopt_keys:
            self.setOption(key,koptions[key])
        #end
        
        
    def __solve__(self, opt_problem={}, *args, **kwargs):
        
        '''
        Run Optimizer (Optimizer Specific Routine)
        
        **Keyword arguments:**
        
        - opt_problem -> INST: Optimization problem instance, *Default* = {}
        
        Documentation last updated:  Feb. 03, 2011 - Peter W. Jansen
        '''
        
        pass
        
        
    def __call__(self, opt_problem={}, *args, **kwargs):
        
        '''
        Run Optimizer (Calling Routine)
        
        **Keyword arguments:**
        
        - opt_problem -> INST: Optimization problem instance, *Default* = {}
        
        Additional arguments and keyword arguments are passed to the objective function call
        
        Documentation last updated:  Feb. 03, 2011 - Peter W. Jansen
        '''
        
        # Check Optimization Problem
        if not isinstance(opt_problem,Optimization):
            try:
                hasattr(opt_problem,'_constraints')
            except:
                raise ValueError("Input is not a Valid Optimization Problem Instance\n")
            #end
        #end
        
        # Solve Optimization Problem
        return self.__solve__(opt_problem, *args, **kwargs)
        
    def _on_setOption(self, name, value):
        
        '''
        Set Optimizer Option Value (Optimizer Specific Routine)
        
        **Arguments:**
        
        - name -> STR: Option name
        - value ->   : Option value
        
        Documentation last updated:  Feb. 07, 2011 - Peter W. Jansen
        '''
        
        raise NotImplementedError()
        
    def setOption(self, name, value=None):
        
        '''
        Set Optimizer Option Value (Calling Routine)
        
        **Arguments:**
        
        - name -> STR: Option Name
        
        **Keyword arguments:**
        
        - value -> FLOAT/INT/BOOL: Option Value, *Default* = None
        
        Documentation last updated:  Feb. 07, 2011 - Peter W. Jansen
        '''
        
        # 
        def_options = self.options['defaults']
        if def_options.has_key(name):
            if (type(value) == def_options[name][0]):
                self.options[name] = [type(value),value]
            else:
                raise IOError('Incorrect ' + repr(name) + ' value type')
            #end
        else:
            raise IOError(repr(name) + ' is not a valid option name')
        #end
        
        # 
        self._on_setOption(name, value)
        
    def _on_getOption(self, name):
        
        '''
        Get Optimizer Option Value (Optimizer Specific Routine)
        
        **Arguments:**
        
        - name -> STR: Option name
        
        Documentation last updated:  Feb. 07, 2011 - Peter W. Jansen
        '''
        
        raise NotImplementedError()
        
    def getOption(self, name):
        
        '''
        Get Optimizer Option Value (Calling Routine)
        
        **Arguments:**
        
        - name -> STR: Option name
        
        Documentation last updated:  Feb. 07, 2011 - Peter W. Jansen
        '''
        
        # 
        def_options = self.options['defaults']
        if def_options.has_key(name):
            return self.options[name][1]
        else:	
            raise IOError(repr(name) + ' is not a valid option name')
        #end
        
        # 
        self._on_getOption(name)
        
    def _on_getInform(self, info):
        
        '''
        Get Optimizer Result Information (Optimizer Specific Routine)
        
        **Arguments:**
        
        - info -> STR: Information key
        
        Documentation last updated:  Feb. 07, 2011 - Peter W. Jansen
        '''
        
        raise NotImplementedError()
        
    def getInform(self, infocode=None):
        
        '''
        Get Optimizer Result Information (Calling Routine)
        
        **Keyword arguments:**
        
        - infocode -> INT: information code key
        
        Documentation last updated:  Feb. 07, 2011 - Peter W. Jansen
        '''

        if (infocode == None):
            return self.informs
        else:
            return self._on_getInform(infocode)
        # end if
        
    def _on_flushFiles(self):
        
        '''
        Flush Output Files (Optimizer Specific Routine)
        
        Documentation last updated:  August. 09, 2009 - Ruben E. Perez
        '''
        
        raise NotImplementedError()
        
    def flushFiles(self):
        
        '''
        Flush Output Files (Calling Routine)
        
        Documentation last updated:  August. 09, 2009 - Ruben E. Perez
        '''
        
        self._on_flushFiles()       
  
        



#==============================================================================
# Optimizer Test
#==============================================================================
if __name__ == '__main__':
    
    # Test Optimizer
    print 'Testing Optimizer...'
    opt = Optimizer()
    opt.ListAttributes()
    
