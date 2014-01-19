#!/usr/bin/env python
'''
pyOpt_variable

Holds the representation of a single pyOptSparse variable 

Copyright (c) 2008-2013 by pyOpt Developers
All rights reserved.

Developers
----------
    Dr. Gaetan K.W. Kenway (GKK)

History
-------
    v. 1.0  - Initial Class Creation (GKK, 2013)
'''
# =============================================================================
# Variable Class
# =============================================================================
class Variable(object):
    '''
    Variable Class Initialization
    '''
    def __init__(self, name, type, value, lower, upper, scale, scalar=False, choices=None):
        
        self.name = name
        self.type = type
        self.scalar = scalar
        if self.type == 'c':
            self.value = value/scale
            self.lower = lower/scale
            self.upper = upper/scale
            self.scale = scale
        elif self.type == 'i':
            self.value = int(value)
            self.lower = lower
            self.upper = upper
        elif self.type == 'd':
            assert choices is not None, 'A discrete variable requires to input an array of choices'
            self.choices = choices
            self.value = self.choices[int(value)]
            self.lower = 0
            self.upper = len(self.choices)
        # end if
        
    def __str__(self):
        '''
        Print Structured List of Variable
        '''
        
        if (self.type == 'd'):
            return ('Name    Type       Value       Lower Bound  Upper Bound\n'+'	 '+str(self.name).center(9) +'%5s	%14f %14.2e %12.2e \n' %(self.type, self.choices[int(self.value)], min(self.choices), max(self.choices)))
        else:
            return ('Name    Type       Value       Lower Bound  Upper Bound\n'+'	 '+str(self.name).center(9) +'%5s	%14f %14.2e %12.2e \n' %(self.type, self.value, self.lower, self.upper))
        # endif
  
