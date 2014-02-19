#!/usr/bin/env python

try:
    from pyIPOPT import IPOPT
    __all__ = ['IPOPT']
except:
    __all__ = []
# end try
