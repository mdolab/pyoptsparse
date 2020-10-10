Developer's Guide
=================

This document provides a high-level description of the software, the structure and layout of the code, and how different parts interact.
If you plan to develop pyOptSparse, please read this to gain some basic understanding of the code.

Code layout
-----------
pyOptSparse uses the object-oriented approach, so several classes are defined and used throughout the code:

-  :meth:`Objective <pyoptsparse.pyOpt_objective.Objective>`, :meth:`Constraint <pyoptsparse.pyOpt_constraint.Constraint>`, and :meth:`Variable <pyoptsparse.pyOpt_variable.Variable>`
   are objects created internally in pyOptSparse to store various values during an optimization.
-  :meth:`Gradient <pyoptsparse.pyOpt_gradient.Gradient>` is used to compute the gradient, either with the various finite difference-type approximations
   or with functions provided by the user.
-  :meth:`Optimizer <pyoptsparse.pyOpt_optimizer.Optimizer>` is the class that all specific optimizer classes inherit from.
-  :meth:`Optimization <pyoptsparse.pyOpt_optimization.Optimization>` stores the optimization problem, i.e. the ``optProb`` object created in the run scripts.
   -  :meth:`Solution <pyoptsparse.pyOpt_solution.Solution>` stores the optimization solution, and is typically the object returned by the optimization call.
-  :meth:`History <pyoptsparse.pyOpt_history.History>` is effectively a wrapper on SqliteDict to store the iteration history of an optimization.

There are also some miscellaneous classes and functions:

-  :meth:`Utils <pyoptsparse.pyOpt_utils>` module which contains utility functions for working with sparse matrices in the pyOptSparse format.
-  ``Error`` and ``pyOptSparseWarning`` classes for formatting error and warning messages in a specific way.
-  ``MPI`` class which will either use the ``mpi4py`` package or its own fake communicator when ``mpi4py`` is not available.
   This allows for cleaner MPI code in the rest of the package, without worrying about whether MPI is available.

Below we will go into a bit more detail on how certain things are done in pyOptSparse.

Optimization problem definition
-------------------------------
This is the part of the code that stores 


Function and gradient computation
---------------------------------

Optimization
------------



Creating new optimizers
-----------------------


Sparse matrices
---------------

MPI handling
------------