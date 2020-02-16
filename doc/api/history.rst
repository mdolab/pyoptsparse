.. _history:

History
-------

Suppose we have an optimization problem with one DV group ``xvars``, one constraint ``con``, and the objective is called ``obj``.
In this case, the history file would have the following layout::

    ├── varInfo
    │   └── xvars
    │       ├── lower
    │       ├── upper
    │       └── scale
    ├── conInfo
    │   └── con
    │       ├── lower
    │       ├── upper
    │       └── scale
    ├── version
    ├── optimizer
    ├── 0
    │   ├── xuser
    │   │   └── xvars
    │   ├── funcs
    │   │   ├── obj
    │   │   └── con
    │   ├── fail
    │   └── isMajor
    ├── 1
    │   ├── xuser
    │   │   └── xvars
    │   ├── funcsSens
    │   │   ├── obj
    │   │   │   └── xvars
    │   │   └── con
    │   │       └── xvars
    │   ├── fail
    │   └── isMajor
    ├── last
    ├── xs
    └── hs

The main optimization history is indexed via call counters, in this example ``0`` and ``1``.
Note that they do not match the major/minor iterations of a given optimizer, since gradient evaluations are stored separate from the function evaluation.

For SNOPT, a number of other values can be requested and stored in each major iteration, such as the feasibility and optimality from the SNOPT print out file.


API
===
.. currentmodule:: pyoptsparse.pyOpt_history

.. autoclass:: History
   :members:
