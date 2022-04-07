.. _slsqp:

SLSQP
=====
SLSQP optimizer is a sequential least squares programming algorithm
which uses the Han-Powell quasi-Newton method with a BFGS update of
the B-matrix and an L1-test function in the step-length algorithm. The
optimizer uses a slightly modified version of Lawson and Hanson's NNLS
nonlinear least-squares solver.

The version provided is the original source code from 1991 by Dieter Kraft.

Options
-------
.. optionstable:: pyoptsparse.pySLSQP.pySLSQP.SLSQP
   :filename: SLSQP_options.yaml


Informs
-------
.. optionstable:: pyoptsparse.pySLSQP.pySLSQP.SLSQP
   :type: informs

API
---
.. currentmodule:: pyoptsparse.pySLSQP.pySLSQP

.. autoclass:: SLSQP
   :members: __call__

