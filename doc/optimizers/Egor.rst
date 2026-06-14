.. _egor:

Egor
====

Egor is a surrogate-based Efficient Global Optimization (EGO) algorithm provided by the
`egobox <https://github.com/relf/EGObox>`_ package which is intalled with:

.. prompt:: bash

   pip install egobox

Egor uses `bayesian optimization <https://en.wikipedia.org/wiki/Bayesian_optimization>`_ techniques 
well-suited to find the global optimum of an expansive-to-evaluate black-box function. 
Basically, it uses a surrogate model to approximate the objective function and 
an infill criterion (aka acquisition function) to guide the search for the optimum.

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
