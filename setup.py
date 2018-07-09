from __future__ import print_function

import os
import sys

# Check if we have numpy:
try:
    from numpy.distutils.misc_util import Configuration
    import numpy.distutils.core
    from numpy.distutils.core import setup
except:
    raise ImportError('pyOptSparse requires numpy version 1.0 or later')

# HACK to make bdist_wheel command usable when using numpy.distutils.core.setup
try:
    from wheel import bdist_wheel
except ImportError:
    if 'bdist_wheel' in sys.argv:
        print("\nThe bdist_wheel option requires the 'wheel' package to be installed.\n"
              "Install it using 'pip install wheel'.")
        sys.exit(-1)
else:
    numpy.distutils.core.numpy_cmdclass['bdist_wheel'] = bdist_wheel.bdist_wheel


if len(sys.argv) == 1:
    print("\nTo install, run: python setup.py install --user\n\n"
          "To build, run: python setup.py build_ext --inplace\n\n"
          "For help on C-compiler options run: python setup.py build --help-compiler\n\n"
          "For help on Fortran-compiler options run: python setup.py build --help-fcompiler\n\n"
          "To specify a Fortran compiler to use run: python setup.py install --user --fcompiler=<fcompiler name>\n\n"
          "For further help run: python setup.py build --help"
      )
    sys.exit(-1)


def configuration(parent_package='', top_path=None):
    config = Configuration(None, parent_package, top_path)
    config.set_options(ignore_setup_xxx_py=True,
                       assume_default_configuration=True,
                       delegate_options_to_subpackages=True,
                       quiet=True)
    config.add_subpackage('pyoptsparse')
    return config


if __name__ == '__main__':

    setup(
        name             = 'pyoptsparse',
        version          = '1.0.0',
        author           = 'Dr. Gaetan Kenway',
        author_email     = 'gaetank@gmail.com',
        maintainer       = 'Dr. Gaetan Kenway',
        maintainer_email = 'gaetank@gmail.com',
        description      = 'Python package for formulating and solving nonlinear constrained optimization problems',
        long_description = 'pyOptSparse is a Python package for formulating and solving nonlinear constrained optimization problems',
        keywords         = 'optimization',
        license          = 'GNU LGPL',
        platforms        = ['Windows','Linux','Solaris','Mac OS-X','Unix'],
        classifiers      = ['Development Status :: 5 - Production/Stable',
                            'Environment :: Console',
                            'Intended Audience :: Science/Research',
                            'Intended Audience :: Developers',
                            'Intended Audience :: Education',
                            'License :: LGPL',
                            'Operating System :: Microsoft :: Windows',
                            'Operating System :: POSIX :: Linux',
                            'Operating System :: Unix',
                            'Operating System :: MacOS',
                            'Programming Language :: Python',
                            'Topic :: Scientific/Engineering',
                            'Topic :: Software Development',
                            'Topic :: Education'],
        configuration = configuration,
    )
