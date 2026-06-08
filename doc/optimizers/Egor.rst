.. _egor:

Egor
====

Egor is a surrogate-based Efficient Global Optimization (EGO) algorithm provided by the
`egobox <https://github.com/rel/EGObox>`_ package.

The pyOptSparse wrapper is derivative-free and targets single-objective, bounded,
continuous design spaces. Constraint values are passed to Egor with the pyOptSparse
constraint convention transformed to :math:`c(x) \le 0`.

Options
-------
.. optionstable:: pyoptsparse.pyEgor.pyEgor.Egor
   :filename: Egor_options.yaml

API
---
.. currentmodule:: pyoptsparse.pyEgor.pyEgor

.. autoclass:: Egor
   :members: __call__
