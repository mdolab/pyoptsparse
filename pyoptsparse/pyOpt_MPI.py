#!/usr/bin/env python
'''
pyOptSparse_MPI

A simple wrapper to MPI that enables pyOptSparse to work without
mpi4py. Only the method from the COMM object that are actually used in
pyOptSparse are included here. 

Copyright (c) 2008-2013 by Dr. Gaetan Kenway
All rights reserved.

Developers:
-----------
- Dr. Gaetan K.W. Kenway (GKK)

History
-------
    v. 1.0  - Initial Class Creation (GKK, 2013)
'''
from __future__ import print_function
import warnings
class COMM(object):
    def __init__(self):
        self.rank = 0
        self.size = 1
    def bcast(self, obj=None, root=0):
        return obj
    def Reduce(self, sendbuf, recvbuf, op, root=0):
        recvbuf = sendbuf.copy()
        return recvbuf
    def allreduce(self, sendobj=None, recvobj=None, op=None):
        return sendobj
    def gather(self, sendobj, recvobj=None, root=0):
        return [sendobj]
    def recv(self, obj=None, source=0, tag=0, status=None):
        return obj
    
class myMPI(object):
    def __init__(self):
        self.COMM_WORLD = COMM()
        self.SUM = 'SUM'
        self.LOR = 'OR'
try:
    from mpi4py import MPI
except:
    warn = 'mpi4py could not be imported. mpi4py is required to use\
 the parallel gradient analysis and parallel objective analysis for\
 non-gradient based optimizers. Continuing using a dummy MPI module\
 from pyOptSparse.'
    warnings.warn(warn)
    MPI = myMPI()
