#!/usr/bin/env python
"""
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
"""
from __future__ import print_function
# =============================================================================
# External Python modules
# =============================================================================
import os
import shelve
import numpy
eps = numpy.finfo(1.0).eps
from .pyOpt_error import Error, pyOptSparseWarning
# =============================================================================
# History Class
# =============================================================================
class History(object):
    """
    Optimizer History Class Initialization. This is essentially a
    thin wrapper around a shelve dictionary to facilitate
    operations with pyOptSparse

    Parameters
    ----------
    fileName : str
       File name for history file

    temp : bool
       Flag to signify that the file should be deleted after it is
       closed

    flag : str
       String of flags to be passed to shelve.open. The only
       useful flag is 'r' which will cause the database to be
       opened in read-only mode. This is often necessary when the
       history file needs to be read from a read-only partition
       during a HPC run job. 
       """
    def __init__(self, fileName, temp=False, flag=''):

        if flag == '':
            # If we are writing, we expliclty remove the file to
            # prevent old keys from "polluting" the new histrory
            if os.path.exists(fileName):
                os.remove(fileName)
            self.db = shelve.open(fileName, protocol=2, writeback=True)
        else:
            if os.path.exists(fileName):
                self.db = shelve.open(fileName, protocol=2, flag=flag)
            else:
                raise Error('The requested history file %s to open in read-only\
                mode does not exist.'% fileName)

        self.temp = temp
        self.fileName = fileName

        # Load any keys it happens to have:
        self.keys = list(self.db.keys())

    def close(self):
        """Close the underlying database"""
        self.db.close()
        if self.temp:
            os.remove(self.fileName)

    def write(self, callCounter, data):
        """This is the main to write data. Basically, we just pass in
        the callCounter, the integer forming the key, and a dictionary
        which will be written to the key"""
        
        # String key to database on disk
        key = '%d'% callCounter

        self.db[key] = data
        self.db['last'] = key
        self.db.sync()
        
    def writeData(self, key, data):
        """
        Write arbitrary key:data value to db
        """
        self.db[key] = data
        self.db.sync()

    def validPoint(self, callCounter, x):
        """
        Determine if callCounter is in the database AND that
        the x matches the x in that point
        """
        key = '%d'% callCounter

        if key not in self.keys:
            return False
        if len(x) == len(self.db[key]['x']):
            if numpy.linalg.norm(x-self.db[key]['x']) < eps:
                return True
            else:
                return False
        else:
            pyOptSparseWarning("The number of x-variables in the\
            history file '%s' was %d, but the number of variables\
            in the current optimization is %d. Continuing without\
            hot starting."% (self.fileName, len(self.db[key]['x']), len(x)))
            return False

    def read(self, callCounter):
        """
        Read data for index 'callCounter'. Note that this
        point should be verified by calling validPoint() first.
        """
        key = '%d'% callCounter
        return self.db[key]
    
    def readData(self, key):
        """
        Read data for generic key 'key'.
        """
        try:
            return self.db[key]
        except KeyError:
            return None
    
    def __del__(self):
        try:
            self.db.close()
            if self.temp:
                os.remove(self.fileName)
        except:
            pass

#==============================================================================
# Optimizer History Test
#==============================================================================
if __name__ == '__main__':
    
    # Test Optimizer History
    print('Testing Optimizer History...')
    
