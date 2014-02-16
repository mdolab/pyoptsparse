#!/usr/bin/env python
"""
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
"""
from .pyOpt_error import Error
INFINITY=1e20
# =============================================================================
# Variable Class
# =============================================================================
class Variable(object):
    """
    Variable Class Initialization
    """
    def __init__(self, name, type, value, lower, upper, scale,
                 scalar=False, choices=None):

        self.name = name
        self.type = type
        self.scalar = scalar
        if self.type == 'c':
            if lower is None:
                self.lower = -INFINITY
            else:
                self.lower = lower*scale

            if upper is None:
                self.upper = INFINITY
            else:
                self.upper = upper*scale

            self.value = value*scale
            self.scale = scale
        elif self.type == 'i':
            self.value = int(value)
            self.lower = lower
            self.upper = upper
        elif self.type == 'd':
            if choices is not None:
                raise Error('A discrete variable requires \
                to input an array of choices')
            self.choices = choices
            self.value = self.choices[int(value)]
            self.lower = 0
            self.upper = len(self.choices)
        
    def __str__(self):
        """
        Print Structured List of Variable
        """

        res = 'Name     Type       Value       '
        res += 'Lower Bound  Upper Bound\n'

        if self.type == 'd':
            res += '	 '
            res += str(self.name).center(15)
            res += '%25s%20f %14.2e %12.2e \n'% (
                self.type, self.choices[int(self.value)],
                min(self.choices), max(self.choices))
        else:
            lower = self.lower
            upper = self.upper
            if self.lower is None:
                lower = -1e20
            if self.upper is None:
                upper = 1e20
            
            res += '	 '
            res += str(self.name).center(9)
            res += '%5s	%14f %14.2e %12.2e \n'% (
                self.type, self.value, lower, upper)
            
        return res
