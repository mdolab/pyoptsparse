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
SNOPT is available for purchase `here 
<http://www.sbsi-sol-optimize.com/asp/sol_snopt.htm>`_. Upon purchase, you should receive a zip file. Within the zip file, there is a folder called ``src``. To use SNOPT with pyoptsparse, paste all files from ``src`` except snopth.f into ``pyoptsparse/pySNOPT/source``.

From v2.0 onwards, only SNOPT v7.7.x is officially supported.
To use pyOptSparse with previous versions of SNOPT, please checkout release v1.2.
The current version of SNOPT being tested is v7.7.5, although v7.7.1 is also expected to work.


Options
-------
Please refer to the SNOPT user manual for a complete listing of options and their default values.
The following are the options set in Python for the wrapper.

- The SNOPT option ``Proximal iterations limit`` has its default value changed to 10000, in order to fully solve the proximal point problem to optimality
- The option ``Save major iteration variables`` is unique to the Python wrapper, and takes a list of values which can be saved at each iteration to the History file.
  Possible values are

  - ``step``
  - ``merit``
  - ``feasibility``
  - ``optimality``
  - ``penalty``
  - ``Hessian``
  - ``slack``
  - ``lambda``
  - ``condZHZ``
- The option ``Start`` corresponds to the value directly passed to the SNOPT kernel, and will be overwritten if another option, e.g. ``Cold start`` is supplied.
- The default values for ``Total integer workspace`` and ``Total real workspace`` are computed internally based on the number of design variables and constraints.
- If SNOPT determines that ``Total character workspace``, ``Total integer workspace``, or ``Total real workspace`` are too small, the Python wrapper will overwrite these values with estimates for the required workspace lengths from SNOPT and initialize the optimizer for a second time.
  SNOPT might still exit with ``80``, but this should automate the storage allocation for most cases.

.. optimizertable:: pyoptsparse.pySNOPT.pySNOPT.SNOPT
   :type: options


Informs
-------
.. optimizertable:: pyoptsparse.pySNOPT.pySNOPT.SNOPT
   :type: informs

API
---
.. currentmodule:: pyoptsparse.pySNOPT.pySNOPT

.. autoclass:: SNOPT
   :members: __call__

