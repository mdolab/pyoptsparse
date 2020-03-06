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
import os, copy
import numpy
from .pyOpt_error import Error
from sqlitedict import SqliteDict
from collections import OrderedDict
eps = numpy.finfo(numpy.float64).eps
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
                # we cast the db to OrderedDict so we do not have to
                # manually close the underlying db at the end
                self.db = OrderedDict(SqliteDict(fileName))
            else:
                raise Error("The requested history file %s to open in "
                            "read-only mode does not exist."% fileName)
            self._processDB()
        else:
            raise Error('The flag argument to History must be \'r\' or \'n\'.')
        self.temp = temp
        self.fileName = fileName
        self.flag = flag

    def close(self):
        """
        Close the underlying database.
        This should only be used in write mode. In read mode, we close the db
        during initialization.
        """
        if self.flag == 'n':
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
        return key in self.keys


    def read(self, callCounter):
        """
        Read data for index 'callCounter'. Note that this
        point should be verified by calling pointExists() first.
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
            if numpy.isclose(xuser,x,atol=eps,rtol=eps).all() and 'funcs' in self.db[key].keys():
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
        x_array = numpy.hstack(x_list)
        return x_array
    
    def _processDB(self):
        """
        Pre-processes the DB file and store various values into class attributes.
        These will be used later when calling self.getXX functions.
        """
        # Load any keys it happens to have:
        self.keys = list(self.db.keys())
        # load info
        self.DVInfo = self.readData('varInfo')
        self.conInfo = self.readData('conInfo')
        self.objInfo = self.readData('objInfo')
        # load names
        self.DVNames = list(self.DVInfo.keys())
        self.conNames = list(self.conInfo.keys())
        self.objName = list(self.objInfo.keys())

        # extract list of callCounters from self.keys
        # this just checks if each key contains only digits, then cast into int
        self.callCounters = sorted([int(x) for x in self.keys if x.isdigit()])

        # extract all information stored in the call counters
        self.iterKeys = set()
        for i in self.callCounters:
            val = self.read(i)
            self.iterKeys.update(val.keys())
        self.iterKeys = list(self.iterKeys)

        # metadata
        self.metadata = self.readData('metadata')

    def getIterKeys(self):
        return copy.deepcopy(self.iterKeys)
    
    def getDVNames(self):
        # only do this if we open the file with 'r' flag
        if self.flag != 'r':
            return
        return copy.deepcopy(self.DVNames)

    def getConNames(self):
        # only do this if we open the file with 'r' flag
        if self.flag != 'r':
            return
        return copy.deepcopy(self.conNames)
    
    def getObjName(self):
        # only do this if we open the file with 'r' flag
        if self.flag != 'r':
            return
        if len(self.objName) == 1:
            return copy.deepcopy(self.objName[0])
        else:
            return copy.deepcopy(self.objName)
    
    def getObjInfo(self, key=None):
        # only do this if we open the file with 'r' flag
        if self.flag != 'r':
            return
        if key is not None:
            return copy.deepcopy(self.objInfo[key])
        else:
            if len(self.objName) == 1:
                return copy.deepcopy(self.objInfo[self.objName[0]])
            else:
                return copy.deepcopy(self.objInfo)
    
    def getConInfo(self,key=None):
        # only do this if we open the file with 'r' flag
        if self.flag != 'r':
            return
        if key is not None:
            if isinstance(key,str):
                return copy.deepcopy(self.conInfo[key])
            elif isinstance(key,list):
                d = OrderedDict()
                for k in key:
                    d[k] = self.conInfo[k]
                return d
        else:
            return copy.deepcopy(self.conInfo)
    
    def getDVInfo(self,key=None):
        # only do this if we open the file with 'r' flag
        if self.flag != 'r':
            return
        if key is not None:
            if isinstance(key,str):
                return copy.deepcopy(self.DVInfo[key])
            elif isinstance(key,list):
                d = OrderedDict()
                for k in key:
                    d[k] = self.DVInfo[k]
                return d
        else:
            return copy.deepcopy(self.DVInfo)
        
    def getMetadata(self):
        # only do this if we open the file with 'r' flag
        if self.flag != 'r':
            return
        return copy.deepcopy(self.metadata)
    
    def _scaleValues(self, name, values):
        """
        This function scales the values, where the factor is extracted from the
        Info dictionaries, according to "name"
        """
        if name in self.objName:
            factor = self.objInfo[name]['scale']
        elif name in self.conNames:
            factor = self.conInfo[name]['scale']
        elif name in self.DVNames:
            factor = self.DVInfo[name]['scale']
        return values * factor

    def getCallCounters(self):
        # only do this if we open the file with 'r' flag
        if self.flag != 'r':
            return
        return copy.deepcopy(self.callCounters)


    def getIterValues(self, names=None, callCounters=None, major=True, scaled=False):
        """
        Parses an existing history file and returns a data dictionary used to post-process optimization results.

        Parameters
        ----------
        names : list or str
            the values of interest, can be the name of any DV, objective or constraint,
            or a list of them. If None, all values are returned.
        
        callCounters : list of ints
            a list of callCounters to extract information from.
            If the callCounter is invalid, i.e. outside the range or is a funcsSens evaluation, then it is skipped.
            If None, all callCounters are looped over.

        major : bool
            flag to specify whether to include only major iterations.
        
        scaled : bool
            flag to specify whether to apply scaling for the values. True means
            to return values that are scaled the same way as the actual optimization.

        Note that regardless of the major flag, failed evaluations are not returned.
        """
        # only do this if we open the file with 'r' flag
        if self.flag != 'r':
            return
        
        allNames = self.DVNames + self.conNames + self.objName
        # cast string input into a single list
        if isinstance(names,str):
            names = [names]
        elif names is None:
            names = allNames
        # error if names isn't either a DV, con or obj
        if not set(names).issubset(set(allNames)):
            raise ValueError("The name provided is not one of DVNames, conNames or objNames")
        
        # set up dictionary to return
        data = {}
        # pre-allocate list for each input
        for name in names:
            data[name] = []

        if callCounters is None:
            callCounters = self.callCounters
        
        for i in callCounters:
            if self.pointExists(i):
                val = self.read(i)
                if 'funcs' in val.keys(): # we have function evaluation
                    if ((major and val['isMajor']) or not major) and not val['fail']:
                        for name in names:
                            if name in self.DVNames:
                                data[name].append(val['xuser'][name])
                            else:
                                data[name].append(val['funcs'][name])
                        
        # reshape inputs and outputs into numpy arrays
        for name in names:
            # only stack if there's more than one value
            if len(callCounters) > 1:
                data[name] = numpy.vstack(data[name])
            else:
                data[name] = data[name][0]
            # scale the values if needed
            if scaled:
                data[name] = self._scaleValues(name,data[name])

        return data

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
    
