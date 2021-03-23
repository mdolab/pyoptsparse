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
We currently test v7.7.7 and v7.7.1.


Options
-------
Please refer to the SNOPT user manual for a complete listing of options and their default values.
The following are a list of

- options which have values changed from the defaults within SNOPT
- options unique to pyOptSparse, implemented in the Python wrapper and not found in SNOPT

.. optionstable:: pyoptsparse.pySNOPT.pySNOPT.SNOPT
   :filename: optimizers/SNOPT_options.yaml

Informs
-------
.. optionstable:: pyoptsparse.pySNOPT.pySNOPT.SNOPT
   :type: informs

API
---
.. currentmodule:: pyoptsparse.pySNOPT.pySNOPT

.. autoclass:: SNOPT
   :members: __call__

