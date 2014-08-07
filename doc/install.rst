.. _install:

Installation
============

Requirements
------------
pyOpt has the following dependencies:

* Python 2.7+ (Python 3.2+)
* Numpy 1.0+
* Swig 1.3+
* c/FORTRAN compiler (compatible with f2py)
* mpi4py

Notes:

* In Windows MinGW is recommended if c/FORTRAN compilers are not available
* In Linux, the python header files (python-dev) are also required.
* Compatibility on Windows 64bit has not been tested 

Building
--------

The recommended approcah to using pyOptSparse is to install
inplace. This does not require root access. From the ``pyoptsparse``
directory run::
    
  >>> python setup.py inplace

To use pyOptSparse in this case, the user should add the path of the
diretory containing ``pyoptsparse`` to the user's ``PYTHONPATH``
enviromental variable. For example, if the pyoptsparse directory is
located at::

  /home/user/hg/pyoptsparse

The required line in the .bashrc file would be::

  export PYTHONPATH=$PYTHONPATH:/home/user/hg/

To install the ``pyOptSparse`` package in a folder on the Python search path 
(usually in a python site-packages or dist-packages folder) run:
    
>>> python setup.py install

This will typically require root access and thus the command actually needs to be:

>>> sudo python setup.py install

Notes:
    
    * You may want to uninstall any previous version of pyOpt before installing a new 
      version, as there may be conflicts.
    * Some optimizers are licensed and their sources are not included with this distribution. 
      To use them, please request their sources from the authors as indicated in the optimizer 
      LICENSE files, and place them in their respective source folders before installing the package.
    * In Windows, if MinGW is used make sure to install for it the C, C++, and Fortran compilers and run:
      
      >>> python setup.py install --compiler=mingw32
      
    * Installing to site-packages/ requires root privileges on Linux.
    * By default pyOpt will attempt to use compilers available on the system. To get a list of 
      available compilers and their corresponding flag on a specific system use:
      
      >>> python setup.py compilers

