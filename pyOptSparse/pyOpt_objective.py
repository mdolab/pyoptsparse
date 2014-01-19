#!/usr/bin/env python
'''
pyOpt_objective

Holds the representation of a pyOptSparse objective

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
# Objective Class
# =============================================================================
class Objective(object):
    
    '''
    Optimization Objective Class

    Parameters
    ----------

    name : str
        Name of this objective
        '''
    
    def __init__(self, name):
        self.name = name
        self.value = 0.0
        self.optimum = 0.0
    def __str__(self):
        
        '''
        Structured Print of Objective
        
        Documentation last updated:  April. 30, 2008 - Peter W. Jansen
        '''
        
        return ( '        Name        Value        Optimum\n'+'	 '+str(self.name).center(9) +'%12g  %12g\n' %(self.value,self.optimum))
    
