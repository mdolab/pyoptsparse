.. _snopt:

SNOPT
=====
SNOPT is a sparse nonlinear optimizer that is particularly useful for
solving large-scale constrained problems with smooth objective
functions and constraints. The algorithm consists of a sequential
quadratic programming (SQP) algorithm that uses a smooth augmented
Lagrangian merit function, while making explicit provision for
infeasibility in the original problem and in the quadratic programming
subproblems. The Hessian of the Lagrangian is approximated using the
BFGS quasi-Newton update.

Installation
------------

Building from source
********************
SNOPT is available for purchase `here 
<http://www.sbsi-sol-optimize.com/asp/sol_snopt.htm>`_. Upon purchase, you should receive a zip file. Within the zip file, there is a folder called ``src``. To use SNOPT with pyoptsparse, paste all files from ``src`` except snopth.f into ``pyoptsparse/pySNOPT/source``.

From v2.0 onwards, only SNOPT v7.7.x is officially supported.
To use pyOptSparse with previous versions of SNOPT, please checkout release v1.2.
We currently test v7.7.7 and v7.7.1.

Installation by conda
*********************

.. _snopt_by_conda:

When installing via conda, all pyoptsparse binaries are pre-compiled and installed as part of the package.
However, the `snopt` binding module cannot be included as part of the package due to license restrictions.

If you are installing via conda and would like to use SNOPT, you will need to build the `snopt` binding module on your own, and inform `pyoptsparse` that it should use that library.

Suppose you have built the binding file, producing ``snopt.cpython-310.so``, living in the folder ``~/snopt-bind``.

To use this module, set the environment variable, ``PYOPTSPARSE_IMPORT_SNOPT_FROM``, e.g.:

.. code-block:: bash

    PYOPTSPARSE_IMPORT_SNOPT_FROM=~/snopt-bind/

This will attempt to load the ``snopt`` binding module from ``~/snopt-bind``. If the module cannot be loaded from this path, a warning will be raised at import time, and an error will be raised if attempting to run the SNOPT optimizer.

Options
-------
Please refer to the SNOPT user manual for a complete listing of options and their default values.
The following are a list of

- options which have values changed from the defaults within SNOPT
- options unique to pyOptSparse, implemented in the Python wrapper and not found in SNOPT

.. optionstable:: pyoptsparse.pySNOPT.pySNOPT.SNOPT
   :filename: SNOPT_options.yaml

Informs
-------
.. optionstable:: pyoptsparse.pySNOPT.pySNOPT.SNOPT
   :type: informs

API
---
.. currentmodule:: pyoptsparse.pySNOPT.pySNOPT

.. autoclass:: SNOPT
   :members: __call__

