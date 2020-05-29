.. _install:

Installation Instructions
=========================

Requirements
------------
pyOpt has the following dependencies:

* Python 2.7+ (Python 3.2+)
* c/FORTRAN compiler (compatible with f2py)

Please make sure these are installed and available for use.
In order to use NSGA2 and NOMAD, SWIG (v1.3+) is also required.
If those optimizers are not needed, then you do not need to install SWIG.
Simply comment out the corresponding lines in ``pyoptsparse/pyoptsparse/setup.py`` so that they are not compiled.
Python dependencies are automatically handled by ``pip``.

.. note::
  * In Linux, the python header files (python-dev) are also required.
  * We do not support operating systems other than Linux.
    If you want to run this on macOS or Windows, you are on your own.

Installation
------------
The easiest and recommended way to install pyOptSparse is with ``pip``.
First clone the repository into a location which is not on the ``$PYTHONPATH``, for example ``~/packages``.
Then in the root ``pyoptsparse`` folder type::

  pip install -e .

For those not using ``conda``, a user install is needed::

  pip install -e . --user

Notes:
    
* You may want to uninstall any previous version of pyOpt before installing a new 
  version, as there may be conflicts.
* Some optimizers are licensed and their sources are not included with this distribution. 
  To use them, please request their sources from the authors as indicated in the optimizer 
  LICENSE files, and place them in their respective source folders before installing the package.
  Refer to specific optimizer pages for additional information.
