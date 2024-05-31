.. _nlpqlp:

NLPQLP
======

NLPQLP is a sequential quadratic programming (SQP) method which solves
problems with smooth continuously differentiable objective function
and constraints. The algorithm uses a quadratic approximation of the
Lagrangian function and a linearization of the constraints. To
generate a search direction a quadratic subproblem is formulated and
solved. The line search can be performed with respect to two
alternative merit functions, and the Hessian approximation is updated
by a modified BFGS formula.

NLPQLP is a proprietary software, which can be obtained `here <https://www.schittkowski.de/numericalsoftware.php>`_.
The supported versions are v4.2.2 and v5.0.3, but other versions may work.

Options
-------
.. optionstable:: pyoptsparse.pyNLPQLP.pyNLPQLP.NLPQLP
   :filename: NLPQLP_options.yaml


Informs
-------
.. optionstable:: pyoptsparse.pyNLPQLP.pyNLPQLP.NLPQLP
   :type: informs

API
---
.. currentmodule:: pyoptsparse.pyNLPQLP.pyNLPQLP

.. autoclass:: NLPQLP
   :members: __call__

