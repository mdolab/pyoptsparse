.. _snopt:

SNOPT
-----
SNOPT is a sparse nonlinear optimizer that is particularly useful for
solving large-scale constrained problems with smooth objective
functions and constraints. The algorithm consists of a sequential
quadratic programming (SQP) algorithm that uses a smooth augmented
Lagrangian merit function, while making explicit provision for
infeasibility in the original problem and in the quadratic programming
subproblems. The Hessian of the Lagrangian is approximated using a
limited-memory quasi-Newton method. 

.. currentmodule:: pyoptsparse.pyoptsparse.pySNOPT.pySNOPT

.. autoclass:: SNOPT
   :members: __call__

