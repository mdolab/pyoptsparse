#!/usr/bin/env python
'''
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
'''
# =============================================================================
# External Python modules
# =============================================================================
import numpy

# =============================================================================
# Constraint Class
# =============================================================================
class Constraint(object):
    
    '''
    Optimization Constraint Class
    '''
    
    def __init__(self, name, linear, wrt, jac, partialReturnOk, lower, upper, scale):
        
        '''
        Constraint Class Initialization
        '''

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

    def __str__(self):
        
        '''
        Print Constraint
        
        Documentation last updated:  April. 30, 2008 - Peter W. Jansen
        '''
        res = ''
        for i in xrange(self.ncon):
            res +='	 '+str(self.name).center(9) + \
                '	  i %15.2e <= %8f <= %8.2e\n' %(
                self.lower[i],self.value[i],self.upper[i])
        # end for

        return res
