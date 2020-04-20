.. _install:

Installation
============

Requirements
------------
pyOpt has the following dependencies:

* Python 2.7+ (Python 3.2+)
* Numpy 1.16+
* Scipy 1.3+
* six 1.13
* sqlitedict 1.6.0
* c/FORTRAN compiler (compatible with f2py)

Please make sure these are installed and available for use.
In order to use NSGA2 and NOMAD, SWIG (v1.3+) is also required.
If those optimizers are not needed, then you do not need to install SWIG.
Simply comment out the corresponding lines in ``pyoptsparse/pyoptsparse/setup.py`` so that they are not compiled.

To facilitate the installation of Python dependencies, there is a file called ``requirements.txt`` which can be used together with ``pip``.
Simply type

.. code-block:: bash

  pip install -r requirements.txt

In the future, we hope to make the package pip-installable so that dependencies can be managed more easily.

.. note::
  * In Windows MinGW is recommended if c/FORTRAN compilers are not available
  * In Linux, the python header files (python-dev) are also required.
  * Compatibility on Windows 64bit has not been tested

Building
--------
The easiest and recommended way to install pyOptSparse is with ``pip``.
First clone the repository into a location which is not on the ``$PYTHONPATH``, for example ``~/packages``.
Then in the root ``pyoptsparse`` folder type::

  pip install -e .

For those not using ``conda``, a user install is needed::

  pip install -e . --user

Two other ways of installing ``pyOptSparse`` are possible, both require running ``setup.py`` manually.
The first approach is to install the package inplace, similar to the ``-e`` flag used above.
This means that no source code is actually copied to a ``site-packages`` folder, meaning
modifying the source code would have an immediate effect, without the need to re-install the code.
From the root directory run::

  >>> python setup.py build_ext --inplace

To use pyOptSparse in this case, the user should add the path of the
root directory to the user's ``$PYTHONPATH`` enviromental variable.
For example, if the pyoptsparse directory is located at::

  /home/<user>/packages/pyoptsparse

The required line in the ``.bashrc`` file would be::

  export PYTHONPATH=$PYTHONPATH:/home/<user>/packages/pyoptsparse

To install the ``pyOptSparse`` package in a folder on the Python search path 
(usually in a python site-packages or dist-packages folder) run::
    
  >>> python setup.py install --user

This will install the package to ``~/.local`` which is typically found
automatically by Python. If a system wide install is desired the
command to run would be (requiring root access)::

  >>> sudo python setup.py install

.. warning::
  Remember to delete the ``build`` directory first when re-building from source.

Notes:
    
* You may want to uninstall any previous version of pyOpt before installing a new 
  version, as there may be conflicts.
* Some optimizers are licensed and their sources are not included with this distribution. 
  To use them, please request their sources from the authors as indicated in the optimizer 
  LICENSE files, and place them in their respective source folders before installing the package.
  Refer to specific optimizer pages for additional information.
* In Windows, if MinGW is used make sure to install for it the C, C++, and Fortran compilers and run:
  
  >>> python setup.py install --compiler=mingw32
  
* By default pyOpt will attempt to use compilers available on the system. To get a list of 
  available compilers and their corresponding flag on a specific system use:
  
  >>> python setup.py build --help-fcompiler

* To see a list of all available ``setup.py`` options for building run 
  
  >>> python setup.py build --help

* In macOS, you may need to force an update of gcc using ‘brew uninstall gcc’ followed by a fresh
  installation of gcc using ‘brew install gcc’ as ‘brew upgrade gcc’ can be insufficient if you
  have recently updated your macOS version.