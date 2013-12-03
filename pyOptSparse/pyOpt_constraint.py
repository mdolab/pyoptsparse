#!/usr/bin/env python
'''
pyOptSparse_constraint

Holds the Python sparse constraint class

Copyright (c) 2013-2013 by Dr. Gaetan Kenway
All rights reserved.
Revision: 1.0   $Date: 19/09/2013 21:00$


Developers:
-----------
- Dr. Gaetan K.W. Kenway (GKK)

History
-------
    v. 1.0  - Initial Class Creation (GKK, 2013)

'''

# =============================================================================
# Standard Python modules
# =============================================================================
import os, sys

# =============================================================================
# External Python modules
# =============================================================================
import numpy
from scipy import sparse
# =============================================================================
# Extension modules
# =============================================================================

# =============================================================================
# Misc Definitions
# =============================================================================
inf = 10.E+20  # define a value for infinity


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
        
        **Arguments:**
         
        See addConGroup in pyOpt_optimization
        
        Documentation last updated:  Feb. 03, 2011 - Peter W. Jansen
        '''
        
        # 
        self.name = name
        self.linear = linear
        self.type = 'i'
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


#==============================================================================
# Constraint Test
#==============================================================================
if __name__ == '__main__':
    
    print 'Testing ...'
    
 
