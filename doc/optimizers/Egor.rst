.. _egor:

Egor
====

Egor is a surrogate-based Efficient Global Optimization (EGO) algorithm provided by the
`egobox <https://github.com/relf/EGObox>`_.

Egor uses `bayesian optimization <https://en.wikipedia.org/wiki/Bayesian_optimization>`_ techniques
well-suited to find the global optimum of an expansive-to-evaluate black-box function.
Basically, it uses a surrogate model to approximate the objective function and
an infill criterion (aka acquisition function) to guide the search for the optimum.

The pyOptSparse wrapper is derivative-free and targets single-objective, bounded,
continuous design spaces. Constraint values are passed to Egor with the pyOptSparse
constraint convention transformed to :math:`c(x) \le 0`.

Installation
------------

Egor is made available through the `egobox <https://pypi.org/project/egobox/>`_ Python package.

.. prompt:: bash

   pip install egobox


Options
-------

Please refer to the Egor help for a complete listing of options and their default values.

.. prompt:: bash

   python
   >>> import egobox as egx
   >>> help(egx.Egor)
   >>> help(egx.GpConfig)

pyoptSparse expects pickable objects while native Egor structures as GpConfig are not pickable.
To workaround this constraint, the pyOptSparse Egor wrapper uses dictionaries which are accepted by EGor
to update the default field values of Egor structures.
Names and default values of the fields are provided in the descriptions below.

.. optionstable:: pyoptsparse.pyEgor.pyEgor.Egor
   :filename: Egor_options.yaml

Informs
-------
.. optionstable:: pyoptsparse.pyEgor.pyEgor.Egor
   :type: informs

API
---
.. currentmodule:: pyoptsparse.pyEgor.pyEgor

.. autoclass:: Egor
   :members: __call__
