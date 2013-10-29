#!/usr/bin/env python
'''
pyOpt_variable

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
'''

__version__ = '$Revision: $'

'''
To Do:
    - variable checking xl < xu
'''

# =============================================================================
# Standard Python modules
# =============================================================================
import os, sys
import pdb

# =============================================================================
# External Python modules
# =============================================================================
#import external

# =============================================================================
# Extension modules
# =============================================================================
#import extension

# =============================================================================
# Misc Definitions
# =============================================================================
inf = 10.E+20  # define a value for infinity


# =============================================================================
# Variable Class
# =============================================================================
class Variable(object):
    
    '''
    Optimization Variable Class
    '''
    
    def __init__(self, name, type='c', value=0.0, *args, **kwargs):
        
        '''
        Variable Class Initialization
        
        **Arguments:**
        
        - name -> STR: Variable Name
        
        **Keyword arguments:**
        
        - type -> STR: Variable Type ('c'-continuous, 'i'-integer, 'd'-discrete), *Default* = 'c'
        - value -> INT/FLOAT: Variable Value, *Default* = 0.0
        - lower -> INT/FLOAT: Variable Lower Value
        - upper -> INT/FLOAT: Variable Upper Value
        - choices -> LIST: Variable Choices
        
        Documentation last updated:  Feb. 07, 2011 - Peter W. Jansen
        '''
        
        # 
        self.name = name
        self.type = type[0].lower()
        if (type[0].lower() == 'c'):
            self.value = float(value)
            self.lower = -inf
            self.upper = inf
            for key in kwargs.keys():
                if (key == 'lower'):
                    self.lower = float(kwargs['lower'])
                #end
                if (key == 'upper'):
                    self.upper = float(kwargs['upper'])
                #end
            #end
        elif (type[0].lower() == 'i'):
            self.value = int(value)
            self.lower = []
            self.upper = []
            for key in kwargs.keys():
                if (key == 'lower'):
                    self.lower = int(kwargs['lower'])
                #end
                if (key == 'upper'):
                    self.upper = int(kwargs['upper'])
                #end
            #end
            if self.lower == []:
                raise IOError('An integer variable requires to input a lower bound value')
            #end
            if self.upper == []:
                raise IOError('An integer variable requires to input an upper bound value')
            #end
        elif (type[0].lower() == 'd'):
            for key in kwargs.keys():
                if (key == 'choices'):
                    self.choices = kwargs['choices']
                else:
                    raise IOError('A discrete variable requires to input an array of choices')
                #end
            #end
            try:
                self.value = self.choices[int(value)]
            except:
                raise IOError('A discrete variable requires the value input to be a integer pointer value of the choices array')
            #end
            self.lower = int(0)
            self.upper = int(len(self.choices))
        else:
            raise IOError('Variable type not understood -- use either c(ontinuous), i(nteger) or d(iscrete)')
        #end
        
        
    def ListAttributes(self):
        
        '''
        Print Structured Attributes List
        
        Documentation last updated:  March. 10, 2008 - Ruben E. Perez
        '''
        
        ListAttributes(self)
        
        
    def __str__(self):
        
        '''
        Print Structured List of Variable
        
        Documentation last updated:  April. 30, 2008 - Peter W. Jansen
        '''
        
        if (self.type == 'd'):
            return ('Name    Type       Value       Lower Bound  Upper Bound\n'+'	 '+str(self.name).center(9) +'%5s	%14f %14.2e %12.2e \n' %(self.type, self.choices[int(self.value)], min(self.choices), max(self.choices)))
        else:
            return ('Name    Type       Value       Lower Bound  Upper Bound\n'+'	 '+str(self.name).center(9) +'%5s	%14f %14.2e %12.2e \n' %(self.type, self.value, self.lower, self.upper))
        #end
    


#==============================================================================
# 
#==============================================================================
def ListAttributes(self):
        
        '''
        Print Structured Attributes List
        
        Documentation last updated:  March. 24, 2008 - Ruben E. Perez
        '''
        
        print '\n'
        print 'Attributes List of: ' + repr(self.__dict__['name']) + ' - ' + self.__class__.__name__ + ' Instance\n'
        self_keys = self.__dict__.keys()
        self_keys.sort()
        for key in self_keys:
            if key != 'name':
                print str(key) + ' : ' + repr(self.__dict__[key])
            #end
        #end
        print '\n'
    


#==============================================================================
# Variable Test
#==============================================================================
if __name__ == '__main__':
    
    print 'Testing ...'
    
    # Test Variable
    var = Variable('x')
    var.ListAttributes()
    
