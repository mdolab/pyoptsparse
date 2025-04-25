.. _install:

Installation Instructions
=========================

Conda
-----
Conda packages are available on ``conda-forge`` and can be installed via

.. prompt:: bash

  conda install -c conda-forge pyoptsparse

This would install pyOptSparse with the built-in optimizers, as well as IPOPT.
If you wish to use optimizers not packaged by ``conda``, e.g. SNOPT, then you must either build the package from source or use the installation script below.
If you have the SNOPT precompiled library available, it is possible to dynamically link it to pyOptSparse following the instructions on the :ref:`SNOPT installation page<snopt_by_conda>`.

Using an installation script
----------------------------
You can build and install pyOptsparse using a `Python script <https://github.com/OpenMDAO/build_pyoptsparse/>`_ developed by the OpenMDAO team.
For usage, see the instruction on the README of the repo.

This script is particularly useful for installing :ref:`IPOPT<ipopt>` and its dependencies.
It can also support SNOPT installation if you have access to the SNOPT source code.

Building from source
--------------------
Requirements
~~~~~~~~~~~~
pyOptSparse has the following dependencies:

* Python 3.7 or 3.8, though other Python 3 versions will likely work
* C and Fortran compilers.
  We recommend ``gcc`` and ``gfortran`` which can be installed via the package manager for your operating system.

Please make sure these are installed and available for use.
Python dependencies are automatically handled by ``pip``, so they do not need to be installed separately.

.. note::
  * In Linux, the python header files (``python-dev``) are also required.
  * **We do not support operating systems other than Linux.**
    For macOS users, the conda package may work out of the box if you do not need any non-default optimizers.
    Also, the installation script by OpenMDAO likely works on macOS.
    For Windows users, a conda package is on the way, if it's not already in the repos.
    This comes with the same disclaimer as the macOS conda package.
    Alternatively, follow the :ref:`conda build instructions<conda build instruction>` below as this will work on any platform.

Installation
~~~~~~~~~~~~
The easiest and recommended way to install pyOptSparse is with ``pip``.
First clone the repository into a location which is not on the ``$PYTHONPATH``, for example ``$HOME/packages/``.
Then in the root ``pyoptsparse`` folder type

.. prompt:: bash

  pip install .

For those not using virtual environments, a user install may be needed

.. prompt:: bash

  pip install . --user

If you plan to modify pyOptSparse, installing with the developer option, i.e. with ``-e``, will save you from re-installing each time you modify the Python code.

.. note::
  Some optimizers are proprietary, and their sources are not distributed with pyOptSparse.
  To use them, please follow the instructions on specific optimizer pages.

To see the list of installed optimizers, use the following:

.. prompt:: bash

  python -c "import pyoptsparse; print(pyoptsparse.list_optimizers())"

Specifying compilers
~~~~~~~~~~~~~~~~~~~~
To specify a non-default compiler (e.g. something other than ``/usr/bin/gcc``), meson recognizes certain `special environment variables <https://mesonbuild.com/Reference-tables.html#compiler-and-linker-selection-variables>`__.
For example, to specify the Intel compilers, simply run

.. prompt:: bash

  FC=$(which ifort) CC=$(which icc) pip install .

.. _install_optview:

Installing OptView
------------------
OptView and OptView-Dash have separate dependencies that must be installed.
To install pyOptSparse including those dependencies, run

.. prompt:: bash

    pip install .[optview]


Testing
-------
pyOptSparse provides a set of unit and regression tests to verify the installation.
To run these tests, first install ``testflo`` which is a testing framework developed by the OpenMDAO team:

.. prompt:: bash

  pip install testflo

Then, in the project root directory, type:

.. prompt:: bash

  testflo . -v

to run all tests.

If there are failed tests, or tests were skipped involving optimizers that should be installed, then refer to the debugging section below.

Debugging Installation Problems
-------------------------------
You may encounter issues such as

.. code-block::

  There was an error importing the compiled slsqp module

The first thing to do is to do a clean install.
This involves the following steps:

#. Uninstall the package via ``pip``
#. If you did a developer install (with ``-e``), check if there are ``.so`` files in the subdirectories, e.g. ``pyoptsparse/pySLSQP``.
   If so, manually delete all ``.so`` files.
#. Remove the ``meson_build`` directory if present.
#. Run ``pip install`` again and test the installation.


If the issue persists, there is probably a linking or runtime issue.
This can be verified by manually importing the compiled library that's causing the issue, for example with:

.. code-block::

  from pyoptsparse.pySLSQP import slsqp


If this throws a ``missing symbol`` error, then there is likely a linking issue at compile time.
If, on the other hand, this throws a ``error while loading shared libraries``, then it's probably a runtime issue with a shared library.
Check to make sure that the ``$LD_LIBRARY_PATH`` is set correctly, for example when running IPOPT.


Update or Uninstall
-------------------
To update pyOptSparse, first delete the ``meson_build`` directory, then update the package using ``git``.
For stability, users are encouraged to stick to tagged releases.
Install the package normally via ``pip``.

To uninstall the package, type

.. prompt:: bash

  pip uninstall pyoptsparse

.. _conda build instruction:

Conda Build Instructions
------------------------
The following instructions explain how to build and install pyOptSparse in a conda environment.
This has the advantage that ``conda`` can be used to install all the necessary dependencies in an isolated and reproducible environment.
Considering how finicky Windows can be with ABI compatibility among various compilers, this is the recommended approach.
The guide will work on any platform, however.

The only build requirement for the build is a working ``conda`` installation as all compilers and dependencies are pulled from the ``conda-forge`` repos, with the exception of a Windows build, which requires Visual Studio 2017 C++ Build Tools.

First, we need to create the ``conda`` environment.
An ``environment.yml`` file is provided in the ``pyoptsparse`` repo:

.. tabs::

  .. code-tab:: bash Linux/OSX

    conda create -y -n pyos-build
    conda activate pyos-build
    conda config --env --add channels conda-forge
    conda config --env --set channel_priority strict

    conda env update -f .github/environment.yml

  .. code-tab:: powershell Windows

    conda create -y -n pyos-build
    conda activate pyos-build
    conda config --env --add channels conda-forge
    conda config --env --set channel_priority strict

    conda env update -f .github\environment.yml

Finally, build the wheel and install it using pip:

.. tabs::

  .. code-tab:: bash Linux/OSX

    # build wheel
    python -m build -n -x .

    # install wheel
    pip install --no-deps --no-index --find-links dist pyoptsparse

  .. code-tab:: powershell Windows

    # set specific compiler flags
    set CC=cl
    set FC=flang
    set CC_LD=link

    # build wheel
    python -m build -n -x .

    # install wheel
    pip install --no-deps --no-index --find-links dist pyoptsparse
