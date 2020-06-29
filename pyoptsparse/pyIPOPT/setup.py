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


def configuration(parent_package="", top_path=None):

    from numpy.distutils.misc_util import Configuration

    config = Configuration("pyIPOPT", parent_package, top_path)

    # Check if we have Ipopt dir....is so assume the user has setup
    # stuff correctly
    add_ipopt = False
    if os.path.exists("pyoptsparse/pyIPOPT/Ipopt"):
        IPOPT_DIR = os.path.join(top_path, "pyoptsparse/pyIPOPT/Ipopt/")
        IPOPT_LIB = os.path.join(top_path, "pyoptsparse/pyIPOPT/Ipopt/lib")
        IPOPT_INC = os.path.join(IPOPT_DIR, "include/coin/")
        add_ipopt = True
    elif os.getenv("IPOPT_INC") is not None and os.getenv("IPOPT_LIB") is not None:
        IPOPT_LIB = os.getenv("IPOPT_LIB")
        IPOPT_INC = os.getenv("IPOPT_INC")
        add_ipopt = True
    elif os.getenv("IPOPT_DIR") is not None:
        IPOPT_DIR = os.getenv("IPOPT_DIR")
        IPOPT_LIB = os.path.join(IPOPT_DIR, "lib")
        IPOPT_INC = os.path.join(IPOPT_DIR, "include/coin-or/")
        add_ipopt = True

    if add_ipopt:
        numpy_include = numpy.get_include()
        FILES = ["src/callback.c", "src/pyipoptcoremodule.c"]
        config.add_extension(
            "pyipoptcore",
            FILES,
            library_dirs=[IPOPT_LIB],
            libraries=["ipopt", "coinmumps", "coinmetis", "dl", "m", "blas", "lapack"],
            extra_link_args=["-Wl,-rpath,%s -L%s" % (IPOPT_LIB, IPOPT_LIB)],
            include_dirs=[numpy_include, IPOPT_INC],
        )
    return config
