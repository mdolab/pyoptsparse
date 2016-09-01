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
    def __init__(self, name, type, value, lower, upper, scale, offset,
                 scalar=False, choices=None):

        self.name = name
        self.type = type
        self.scalar = scalar
        self.choices = None
        if self.type == 'c':
            if lower is None:
                self.lower = -INFINITY
            else:
                self.lower = (lower-offset)*scale

            if upper is None:
                self.upper = INFINITY
            else:
                self.upper = (upper-offset)*scale

            self.value = (value-offset)*scale
            self.scale = scale
            self.offset = offset
        elif self.type == 'i':
            self.value = int(value)
            self.lower = lower
            self.upper = upper
            self.scale = scale
        elif self.type == 'd':
            if choices is None:
                raise Error("A discrete variable requires "
                            "to input an array of choices.")
            self.choices = choices
            self.value = self.choices[int(value)]
            self.lower = 0
            self.upper = len(self.choices)
            self.scale = scale

    def __eq__(self, other):
        """
        Compare two variable objects
        """
        if (self.name == other.name and self.type == other.type and
            self.scalar == other.scalar and self.upper == other.upper and
            self.lower == other.lower and self.choices == other.choices):
            return True
        else:
            return False

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
