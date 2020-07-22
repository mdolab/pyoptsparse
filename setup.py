import os
import sys

import setuptools  # magic import to allow us to use entry_point

# Check if we have numpy:
try:
    from numpy.distutils.misc_util import Configuration
    import numpy.distutils.core
    from numpy.distutils.core import setup
except:
    raise ImportError("pyOptSparse requires numpy version 1.0 or later")

# HACK to make bdist_wheel command usable when using numpy.distutils.core.setup
try:
    from wheel import bdist_wheel
except ImportError:
    if "bdist_wheel" in sys.argv:
        print(
            (
                "\nThe bdist_wheel option requires the 'wheel' package to be installed.\n"
                + "Install it using 'pip install wheel'."
            )
        )
        sys.exit(-1)
else:
    numpy.distutils.core.numpy_cmdclass["bdist_wheel"] = bdist_wheel.bdist_wheel


if len(sys.argv) == 1:
    print(
        (
            "\nTo install, run: python setup.py install --user\n\n"
            + "To build, run: python setup.py build_ext --inplace\n\n"
            + "For help on C-compiler options run: python setup.py build --help-compiler\n\n"
            + "For help on Fortran-compiler options run: python setup.py build --help-fcompiler\n\n"
            + "To specify a Fortran compiler to use run: python setup.py install --user --fcompiler=<fcompiler name>\n\n"
            + "For further help run: python setup.py build --help"
        )
    )
    sys.exit(-1)


def configuration(parent_package="", top_path=None):
    config = Configuration(None, parent_package, top_path)
    config.set_options(
        ignore_setup_xxx_py=True, assume_default_configuration=True, delegate_options_to_subpackages=True, quiet=True
    )
    config.add_subpackage("pyoptsparse")
    return config


if __name__ == "__main__":

    import re

    __version__ = re.findall(r"""__version__ = ["']+([0-9\.]*)["']+""", open("pyoptsparse/__init__.py").read(),)[0]

    setup(
        name="pyoptsparse",
        version=__version__,
        description="Python package for formulating and solving nonlinear constrained optimization problems",
        long_description="pyOptSparse is a Python package for formulating and solving nonlinear constrained optimization problems",
        keywords="optimization",
        install_requires=["sqlitedict>=1.6.0", "numpy>=1.16.4", "scipy>1.2.1", "six>=1.13"],
        platforms=["Linux"],
        classifiers=[
            "Development Status :: 5 - Production/Stable",
            "Environment :: Console",
            "Intended Audience :: Science/Research",
            "Intended Audience :: Developers",
            "Intended Audience :: Education",
            "License :: OSI Approved :: GNU Library or Lesser General Public License (LGPL)",
            "Operating System :: POSIX :: Linux",
            "Programming Language :: Python",
            "Topic :: Scientific/Engineering",
            "Topic :: Software Development",
            "Topic :: Education",
        ],
        configuration=configuration,
        entry_points={"gui_scripts": ["optview = pyoptsparse.postprocessing.OptView:main"]},
    )
