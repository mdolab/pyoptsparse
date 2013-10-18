#!/usr/bin/env python
'''
pyOpt_history

Holds the Python Design Optimization History Class.

Copyright (c) 2008-2013 by pyOpt Developers
All rights reserved.
Revision: 1.0   $Date: 11/12/2009 21:00$


Developers:
-----------
- Dr. Gaetan K. W. Kneway (GKK)

History
-------
    v. 1.0  - Initial Class Creation (GKK 2013)
'''

__version__ = '$Revision: $'

'''
To Do:

'''

# =============================================================================
# Standard Python modules
# =============================================================================
import os, sys
import shelve

# =============================================================================
# External Python modules
# =============================================================================
import numpy
eps = numpy.finfo(1.0).eps
# =============================================================================
# History Class
# =============================================================================
class History(object):
    
    '''
    Abstract Class for Optimizer History Object
    '''
    
    def __init__(self, fileName, temp=False):
        '''
        Optimizer History Class Initialization
        
        **Arguments:**
        - fileName  -> STR: Name for .bin and .cue file
        - temp -> BOOL: If this is a temporary file, delete file
                        upon closing.
        '''
        
        # 
        self.db = shelve.open(fileName, protocol=2, writeback=True)
        self.temp = temp
        self.fileName = fileName

        # Load any keys it happens to have:
        self.keys = self.db.keys()

    def close(self):
        self.db.close()
        if self.temp:
            os.remove(self.fileName)
        # end if

    def write(self, callCounter, fobj, fcon, fail, xn, x_array, gradEvaled, gobj, gcon, **kwargs):
        data = {'fobj':fobj, 'fcon':fcon, 'fail':fail,'x':xn,'x_array':x_array}
        if gradEvaled:
            data['gobj'] = gobj
            data['gcon'] = gcon

        # Add any extra keyword arguments that specific optimizers may
        # want to save
        data.update(kwargs)

        # String key to database on disk
        key = '%d'%callCounter

        self.db[key] = data
        self.db['last'] = key
        self.db.sync()
        
        return

    def validPoint(self, callCounter, x):
        '''
        Determine if callCounter is in the database AND that
        the x matches the x in that point
        '''
        key = '%d'%callCounter     
        if key not in self.keys:
            return False
        # end if
        x_array = self.db[key]['x_array']
        
        # Check that the x's *actually* match
        diff = numpy.dot(x-x_array, x-x_array)

        if diff < eps:
            return True
        else:
            return False
        # end if

    def read(self, callCounter):
        '''
        Read data for index 'callCounter'. Note that this
        point should be verified by calling validPoint() first.'''
        key = '%d'%callCounter
        return self.db[key]

    def __del__(self):
        self.close()

#==============================================================================
# Optimizer History Test
#==============================================================================
if __name__ == '__main__':
    
    # Test Optimizer History
    print 'Testing Optimizer History...'
    
