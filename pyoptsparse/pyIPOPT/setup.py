# Originally contributed by Lorne McIntosh.
# Modified by Eric Xu
# Further modification by random internet people.

# You will probably have to edit this file in unpredictable ways
# if you want pyipopt to work for you, sorry.

# When I installed Ipopt from source, I used the
# --prefix=/usr/local
# option, so this is where I want pyipopt to look for my ipopt installation.
# I only installed from source because the ipopt packaging
# for my linux distribution was buggy,
# so by the time you read this the bugs have probably been fixed
# and you will want to specify a different directory here.

# Futher modification by Gaetan Kenway to work with the pyOptSparse
# build system.

import os, sys
import numpy

def configuration(parent_package='',top_path=None):

    from numpy.distutils.misc_util import Configuration

    numpy_include = numpy.get_include()
    IPOPT_DIR = './Ipopt/'
    IPOPT_LIB = './Ipopt/lib/'
    IPOPT_INC = os.path.join(IPOPT_DIR, 'include/coin/')
    FILES = ['src/callback.c', 'src/pyipoptcoremodule.c']
    config = Configuration('pyIPOPT', parent_package, top_path)
    config.add_extension('pyipoptcore',
                         FILES,
                         library_dirs=[IPOPT_LIB],
                         libraries=['ipopt', 'coinblas','coinhsl','coinlapack','dl','m'],
                         extra_link_args=['-Wl,--rpath -Wl,%s -L%s'% (IPOPT_LIB,IPOPT_LIB)],
                         include_dirs=[numpy_include, IPOPT_INC])
    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
