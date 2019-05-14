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
<http://www.sbsi-sol-optimize.com/asp/sol_snopt.htm>`_. Upon purchase, you should receive a zip file. Within the zip file, there is a folder called "src". To use SNOPT with pyoptsparse, paste all files from `src` except snopth.f into `pyoptsparse/pySNOPT/source`.


API
---


.. currentmodule:: pyoptsparse.pySNOPT.pySNOPT

.. autoclass:: SNOPT
   :members: __call__

