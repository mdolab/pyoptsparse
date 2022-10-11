.. _nsga2:

NSGA2
=====

This optimizer is a non-dominating sorting genetic algorithm that
solves non-convex and non-smooth single and multiobjective
optimization problems. The algorithm attempts to perform global
optimization, while enforcing constraints using a tournament
selection-based strategy

.. warning:: Currently, the Python wrapper does not catch
  exceptions. If there is **any** error in the user-supplied function,
  you will get a seg-fault and no idea where it happened. Please make
  sure the objective is without errors before trying to use nsga2.

Options
-------
.. optionstable:: pyoptsparse.pyNSGA2.pyNSGA2.NSGA2
   :filename: NSGA2_options.yaml


API
---
.. currentmodule:: pyoptsparse.pyNSGA2.pyNSGA2

.. autoclass:: NSGA2
   :members: __call__

