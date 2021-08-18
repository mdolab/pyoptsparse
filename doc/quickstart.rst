.. _quickstart:

Quickstart
==========

The following shows how to get started with pyOptSparse by solving Schittkowski’s TP37 constrained problem.
First, we show the complete program listing and then go through each statement line by line:

.. literalinclude:: ../examples/tp037.py
   :start-after: # rst begin import
   :end-before: # rst end opt

Start by importing the pyOptSparse package:

.. literalinclude:: ../examples/tp037.py
   :start-after: # rst begin import
   :end-before: # rst begin objfunc

Next we define the objective function that takes in the design variable *dictionary* and returns
a *dictionary* containing the constraints and objective, as well as a (boolean) flag indicating if
the objective function evaluation was successful.
For the TP37, the objective function is a simple analytic function:

.. literalinclude:: ../examples/tp037.py
   :start-after: # rst begin objfunc
   :end-before: # rst begin optProb

Notes:

*   The ``xdict`` variable is a dictionary whose keys are the names from each 
    :meth:`addVar <pyoptsparse.pyOpt_optimization.Optimization.addVar>` and 
    :meth:`addVarGroup <pyoptsparse.pyOpt_optimization.Optimization.addVarGroup>` call. The line

    .. code-block:: python

        x = xdict["xvars"]

    retrieves an array of length 3 which are all the variables for this optimization. 

*   The line

    .. code-block:: python

        conval = [0] * 2

    creates a list of length 2, which stores the numerical values of the two constraints.
    The ``funcs`` dictionary return must contain keys that match the constraint names from
    :meth:`addCon <pyoptsparse.pyOpt_optimization.Optimization.addCon>` and
    :meth:`addConGroup <pyoptsparse.pyOpt_optimization.Optimization.addConGroup>`
    as well as the objectives from :meth:`addObj <pyoptsparse.pyOpt_optimization.Optimization.addObj>` calls.
    This is done in the following calls:

    .. code-block:: python

        funcs["obj"] = -x[0] * x[1] * x[2]
        funcs["con"] = conval

Now the optimization problem can be initialized:

.. literalinclude:: ../examples/tp037.py
    :start-after: # rst begin optProb
    :end-before: # rst begin addVar

This creates an instance of the optimization class with a name and a reference to the objective function.
To complete the setup of the optimization problem, the design variables and constraints need to be defined. 

Design variables and constraints can be added either one-by-one or as a group.
Adding variables by group is generally recommended for related variables:

.. literalinclude:: ../examples/tp037.py
    :start-after: # rst begin addVar
    :end-before: # rst begin addCon

This calls adds a group of 3 variables with name ``xvars``.
The variable bounds (side constraints) are 0 for the lower bounds, and 42 for the upper bounds.
The initial values for each variable is 10.0.

Now, we must add the constraints.
Like design variables, these may be added individually or by group.
It is recommended that related constraints are added by group where possible:

.. literalinclude:: ../examples/tp037.py
    :start-after: # rst begin addCon
    :end-before: # rst begin addObj

This call adds two variables with name ``con``.
There is no lower bound for the variables and the upper bound is 0.0. 

We must also assign the the key value for the objective using the
:meth:`addObj <pyoptsparse.pyOpt_optimization.Optimization.addObj>` call:

.. literalinclude:: ../examples/tp037.py
   :start-after: # rst begin addObj
   :end-before: # rst begin print

The optimization problem can be printed to verify that it is set up correctly:

.. literalinclude:: ../examples/tp037.py
   :start-after: # rst begin print
   :end-before: # rst begin OPT

which produces the following table::

   Optimization Problem -- TP037 Constraint Problem
   ================================================================================
      Objective Function: objfunc


      Objectives
         Index  Name            Value
            0  obj     0.000000E+00

      Variables (c - continuous, i - integer, d - discrete)
         Index  Name      Type      Lower Bound            Value      Upper Bound     Status
            0  xvars_0      c     0.000000E+00     1.000000E+01     4.200000E+01           
            1  xvars_1      c     0.000000E+00     1.000000E+01     4.200000E+01           
            2  xvars_2      c     0.000000E+00     1.000000E+01     4.200000E+01           

      Constraints (i - inequality, e - equality)
         Index  Name Type          Lower           Value           Upper    Status  Lagrange Multiplier (N/A)
            0  con    i  -1.000000E+20    0.000000E+00    0.000000E+00         u    9.00000E+100
            1  con    i  -1.000000E+20    0.000000E+00    0.000000E+00         u    9.00000E+100

To solve an optimization problem with pyOptSparse an optimizer must be initialized.
The initialization of one or more optimizers is independent of the initialization of the optimization problem.
To initialize ``SLSQP``, which is an open-source, sequential least squares programming algorithm that comes as part of the
pyOptSparse package, use:

.. literalinclude:: ../examples/tp037.py
   :start-after: # rst begin OPT
   :end-before: # rst begin solve

This initializes an instance of ``SLSQP`` with the option ``IPRINT`` set to -1.
All other options will be set to the default values, which can be found on the optimizer-specific pages.
For example, the default options for SLSQP can be found in the page :ref:`SLSQP`.

Now TP37 can be solved using pyOptSparse’s automatic finite difference for the gradients:

.. literalinclude:: ../examples/tp037.py
   :start-after: # rst begin solve
   :end-before: # rst begin check

We can print the solution objection to view the result of the optimization:

.. literalinclude:: ../examples/tp037.py
   :start-after: # rst begin check
   :end-before: # rst end opt

which produces the following output::

   Optimization Problem -- TP037 Constraint Problem
   ================================================================================
      Objective Function: objfunc

      Solution: 
   --------------------------------------------------------------------------------
      Total Time:                    0.0062
         User Objective Time :       0.0001
         User Sensitivity Time :     0.0011
         Interface Time :            0.0046
         Opt Solver Time:            0.0003
      Calls to Objective Function :      22
      Calls to Sens Function :            8


      Objectives
         Index  Name            Value
            0  obj    -3.456000E+03

      Variables (c - continuous, i - integer, d - discrete)
         Index  Name      Type      Lower Bound            Value      Upper Bound     Status
            0  xvars_0      c     0.000000E+00     2.399997E+01     4.200000E+01           
            1  xvars_1      c     0.000000E+00     1.200001E+01     4.200000E+01           
            2  xvars_2      c     0.000000E+00     1.200000E+01     4.200000E+01           

      Constraints (i - inequality, e - equality)
         Index  Name Type          Lower           Value           Upper    Status  Lagrange Multiplier (N/A)
            0  con    i  -1.000000E+20    7.564591E-07    0.000000E+00         u    9.00000E+100
            1  con    i  -1.000000E+20   -7.200000E+01    0.000000E+00              9.00000E+100
