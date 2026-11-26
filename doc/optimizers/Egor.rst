.. _egor:

Egor
====

Egor is the surrogate-based Efficient Global Optimization (EGO) optimizer from the open-source
`EGObox <https://github.com/relf/EGObox>`_ library.

Egor uses `bayesian optimization <https://en.wikipedia.org/wiki/Bayesian_optimization>`_ techniques
well-suited to find the global optimum of an expansive-to-evaluate black-box function.
Basically, it uses a surrogate model to approximate the objective function and
an infill criterion (aka acquisition function) to guide the search for the optimum.

The pyOptSparse wrapper is derivative-free and targets single-objective, bounded,
continuous design spaces: all design variables require finite lower and upper bounds,
otherwise a ``ValueError`` listing the unbounded ones is raised.
Constraint values are passed to Egor as they are, with their pyOptSparse bounds given as
constraint specs (``leq``, ``geq``, ``eq`` or ``between``) and tolerances set with the ``cstr_tols`` option.

Installation
------------

Egor is made available through the `egobox <https://pypi.org/project/egobox/>`_ Python package.

.. prompt::

   $ pip install egobox

``egobox`` is also available via conda-forge::

   $ conda install -c conda-forge egobox


Options
-------

Please refer to the Egor help for a complete listing of options and their default values.

.. prompt::

   $ python
   >>> import egobox as egx
   >>> help(egx.Egor)
   >>> help(egx.GpConfig)

pyoptSparse expects pickable objects while native Egor structures as GpConfig are not pickable.
To workaround this constraint, the pyOptSparse Egor wrapper uses dictionaries which are accepted by Egor
to update the default field values of Egor structures.
Names and default values of the fields are provided in the descriptions below.

.. optionstable:: pyoptsparse.pyEgor.pyEgor.Egor
   :filename: Egor_options.yaml
