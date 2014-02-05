#!/usr/bin/env python
"""
pyOpt_constraint

Holds the representation of a pyOptSparse constraint group

Copyright (c) 2008-2013 by pyOpt Developers
All rights reserved.

Developers:
-----------
- Dr. Gaetan K.W. Kenway (GKK)

History
-------
    v. 1.0  - Initial Class Creation (GKK, 2013)
"""
# =============================================================================
# External Python modules
# =============================================================================
import numpy
eps = numpy.finfo(1.0).eps
# =============================================================================
# Constraint Class
# =============================================================================
class Constraint(object):
    """
    Constraint Class Initialization
    """
    def __init__(self, name, linear, wrt, jac, partialReturnOk,
                 lower, upper, scale):

        self.name = name
        self.linear = linear
        self.wrt = wrt
        self.jac = jac
        self.partialReturnOk = partialReturnOk
        self.lower = lower
        self.upper = upper
        self.value = numpy.zeros_like(lower)
        self.ncon = len(upper)
        self.rs = None
        self.re = None
        self.scale = scale

        # Now we determine what kind of constraint this is: 
        # 1. An equality constraint
        # 2. A upper bound on a 1-sided constraint
        # 3. A lower bound on a 1-sided constraint
        # 4. Lower and Upper bounds on 2-sided constraint
        
        # The first 3, will give a single "constraint" in all
        # optimizers. Some optimizers can only do 1-sided constraints
        # so type 4 must be split into two separate constraints
        # automatically. 

        # Type 1: Equality



    def __str__(self):
        """
        Print Constraint
        
        Documentation last updated:  April. 30, 2008 - Peter W. Jansen
        """
        res = ''
        for i in range(self.ncon):
            res += '	 '+str(self.name).center(9) + \
                   '	  i %15.2e <= %8f <= %8.2e\n' %(
                self.lower[i],self.value[i],self.upper[i])
       
        return res
