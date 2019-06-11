#!/usr/bin/env python
"""
pyOpt_history

Holds the Python Design Optimization History Class.

Copyright (c) 2008-2013 by pyOpt Developers
All rights reserved.
Revision: 1.0   $Date: 11/12/2009 21:00$

Developers:
-----------
- Dr. Gaetan K. W. Kenway (GKK)

History
-------
    v. 1.0  - Initial Class Creation (GKK 2013)
"""
from __future__ import print_function
# =============================================================================
# External Python modules
# =============================================================================
import os
import numpy as np
from .pyOpt_error import Error
from .sqlitedict.sqlitedict import SqliteDict
# =============================================================================
# History Class
# =============================================================================
class History(object):
    """
    Optimizer History Class Initialization. This is essentially a
    thin wrapper around a SqliteDict dictionary to facilitate
    operations with pyOptSparse

    Parameters
    ----------
    fileName : str
       File name for history file

    temp : bool
       Flag to signify that the file should be deleted after it is
       closed

    flag : str
        String specifying the mode. Similar to what was used in
        shelve. 'n' for a new database and 'r' to read an existing one. 
       """
    def __init__(self, fileName, temp=False, flag='n'):

        if flag == 'n':
            # If writing, we expliclty remove the file to
            # prevent old keys from "polluting" the new histrory
            if os.path.exists(fileName):
                os.remove(fileName)
            self.db = SqliteDict(fileName)
        elif flag == 'r':
            if os.path.exists(fileName):
                self.db = SqliteDict(fileName)
            else:
                raise Error("The requested history file %s to open in "
                            "read-only mode does not exist."% fileName)
        else:
            raise Error('The flag argument to History must be \'r\' or \'n\'.')
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
        # if the point exists, we merely update with new data
        if self.pointExists(callCounter):
            oldData = self.read(callCounter)
            oldData.update(data)
            self.db[key] = oldData
        else:
            self.db[key] = data
        self.db['last'] = key
        self.db.sync()
        self.keys = list(self.db.keys())
        
    def writeData(self, key, data):
        """
        Write arbitrary key:data value to db
        """
        self.db[key] = data
        self.db.commit()
        self.keys = list(self.db.keys())

    def pointExists(self, callCounter):
        """
        Determine if callCounter is in the database
        """
        key = '%d'% callCounter

        if key in self.keys:
            return True
        else:
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
    
    def getCallCounter(self,x):
        """
        Returns the callCounter corresponding to the function evaluation at 'x',
        returns None if the point did not match previous evaluations
        """
        last = int(self.db['last'])
        callCounter = None
        for i in range(last,0,-1):
            key = '%d'% i
            xuser = self.deProcessX(self.db[key]['xuser'])
            if np.isclose(x,xuser).all() and 'funcs' in self.db[key].keys():
                callCounter = i
                break
        return callCounter

    def deProcessX(self,xuser):
        """
        This is a much more simple version of pyOpt_history.deProcessX without error checking.
        We traverse the OrderedDict and stack all the DVs as a single numpy array, preserving 
        the order so that we get the correct x vector.
        """
        x_list = []
        for key in xuser.keys():
            x_list.append(xuser[key])
        x_array = np.hstack(x_list)
        return x_array

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
    
