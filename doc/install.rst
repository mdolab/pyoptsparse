.. _install:

Installation Instructions
=========================

Requirements
------------
pyOptSparse has the following dependencies:

* Python 2.7, 3.7 or 3.8, though other Python 3 versions will likely work
* C/FORTRAN compiler (compatible with f2py)

Please make sure these are installed and available for use.
In order to use NSGA2 and NOMAD, SWIG (v1.3+) is also required.
If those optimizers are not needed, then you do not need to install SWIG.
Simply comment out the corresponding lines in ``pyoptsparse/pyoptsparse/setup.py`` so that they are not compiled.
The corresponding lines in ``pyoptsparse/__init__.py`` must be commented out as well.

Python dependencies are automatically handled by ``pip``, so they do not need to be installed separately.

.. note::
  * In Linux, the python header files (python-dev) are also required.
  * We do not support operating systems other than Linux.
    If you want to run this on macOS or Windows, you are on your own.

Installation
------------
The easiest and recommended way to install pyOptSparse is with ``pip``.
First clone the repository into a location which is not on the ``$PYTHONPATH``, for example ``$HOME/packages/``.
Then in the root ``pyoptsparse`` folder type::

  pip install .

For those not using virtual environments, a user install may be needed::

  pip install . --user

If you plan to modify pyOptSparse, installing with the developer option, i.e. with ``-e``, will save you from re-installing each time you modify the Python code.

It is also possible to install pyOptSparse by calling ``python setup.py install``, but this is not recommended.

.. note::
  Some optimizers are proprietary and their sources are not distributed with pyOptSparse.
  To use them, please follow the instructions on specific optimizer pages.

Update or Uninstall
-------------------
To update pyOptSparse, first delete the ``build`` directory, then update the package using ``git``.
For stability, users are encouraged to stick to tagged releases.
Install the package normally via ``pip``.

To uninstall the package, type::

  pip uninstall pyoptsparse

.. note::
  pyOptSparse can optionally run in parallel if a suitable ``mpi4py``
  installation exists. This will be automatically detected and
  imported at run-time.

  If you only want to run in parallel, you can
  force pyOptSparse to do so by setting the environment variable
  ``PYOPTSPARSE_REQUIRE_MPI`` to anyone of these values: ``['always', '1', 'true', 'yes']``
  If a suitable ``mpi4py`` is not available, an exception will be raised and the run
  terminated.

  If you explicitly do not wish to use ``mpi4py``, set the environment variable ``PYOPTSPARSE_REQUIRE_MPI``
  to anything other than those values. This can come in handy, for example, if your ``MPI`` installation
  is not functioning properly, but you still need to run serial code.

