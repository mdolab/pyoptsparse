.. _nlpql:

NLPQL
-----

NLPQL is a sequential quadratic programming (SQP) method which solves
problems with smooth continuously differentiable objective function
and constraints. The algorithm uses a quadratic approximation of the
Lagrangian function and a linearization of the constraints. To
generate a search direction a quadratic subproblem is formulated and
solved. The line search can be performed with respect to two
alternative merit functions, and the Hessian approximation is updated
by a modified BFGS formula. 

.. currentmodule:: pyoptsparse.pyoptsparse.pyNLPQL.pyNLPQL

.. autoclass:: NLPQL
   :members: __call__

