#!/usr/bin/env python
'''
pyOptSparse_error

Holds a simple error handling class for pyOptSparse

Copyright (c) 2008-2013 by Dr. Gaetan Kenway
All rights reserved.

Developers:
-----------
- Dr. Gaetan K.W. Kenway (GKK)

History
-------
    v. 1.0  - Initial Class Creation (GKK, 2013)
'''

class Error(Exception):
    """General error raised when we detect something that the user
    should not have done"""

    def __init__(self, message):
        print '+'+'-'*78+'+'
        print '| pyOptSparse Error: ',
        i = 22
        aux = message.split()
        for word in aux:
            if len(word) + i > 78:
                print ' '*(79-i)+'|'
                print '|',
                i = 2
                print word,
                i += len(word)+1
            else:
                print word,
                i += len(word)+1
            # end if
        # end for
        print ' '*(79-i)+'|'
        print '+'+'-'*78+'+'
