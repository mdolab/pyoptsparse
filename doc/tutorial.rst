.. _tutorial:

Tutorial
========

The following shows how to get started with pyOptSparse by solving
Schittkowski’s TP37 constrained problem. First, we show the complete
program listing and then go through each statement line by line::

  import pyoptsparse
  def objfunc(xdict):
      x = xdict['xvars']
      funcs = {}
      funcs['obj'] = -x[0]*x[1]*x[2]
      conval = [0]*2
      conval[0] = x[0] + 2.*x[1] + 2.*x[2] - 72.0
      conval[1] = -x[0] - 2.*x[1] - 2.*x[2]
      funcs['con'] = conval
      fail = False

      return funcs, fail

  optProb = pyoptsparse.Optimization('TP037', objfunc)
  optProb.addVarGroup('xvars',3, 'c',lower=[0,0,0], upper=[42,42,42], value=10)
  optProb.addConGroup('con',2, lower=None, upper=0.0)
  optProb.addObj('obj')
  print optProb
  opt = pyoptsparse.SLSQP()
  sol = opt(optProb, sens='FD')
  print sol

Start by importing the pyOptSparse package::

  >>> import pyoptsparse

Next we define the objective function that takes in the design
variable *dictionary* and returns a *dictionary* containing the
constraints and objective, as well as a (boolean) flag indicating if
the objective function evaluation was successful. For the TP37, the
objective function is a simple analytic function::

  def objfunc(xdict):
      x = xdict['xvars']
      funcs = {}
      funcs['obj'] = -x[0]*x[1]*x[2]
      conval = [0]*2
      conval[0] = x[0] + 2.*x[1] + 2.*x[2] - 72.0
      conval[1] = -x[0] - 2.*x[1] - 2.*x[2]
      funcs['con'] = conval
      fail = False

      return funcs, fail

Notes:

  1. The ``xdict`` variable is a dictionary whose keys are the names
     from each ``addVar()`` and ``addVarGroup()`` call. The line::

       x = xdict['xvars']

     retrieves an array of length 3 which are all the variables for
     this optimization. 

  2. The line::
    
       conval = [0]*2

     creates a list of length 2, which stores the numerical values of
     the two constraints. The ``funcs`` dictionary return must contain
     keys that match the constraint names from ``addCon`` and
     ``addConGroup``  as well as the objectives from ``addObj`` calls. This 
     is done in the following calls::

       funcs['obj'] = -x[0]*x[1]*x[2]
       funcs['con'] = conval

Now the optimization problem can be initialized::

  >>> optProb = Optimization('TP037', objfunc)

This creates an instance of the optimization class with a name and a
reference to the objective function. To complete the setup of the
optimization problem, the design variables and constraints need to be defined. 

Design variables and constraints can be added either one-by-one or as
a group. Adding variables by group is generally recommended for
related variables::

  >>> optProb.addVarGroup('xvars', 3, 'c', lower=[0,0,0], upper=[42,42,42], value=10)

This calls adds a group of 3 variables with name 'xvars'. The variable
bounds (side constraints) are 0 for the lower bounds, 42 for the upper
bounds. The inital values for each variable is 10.0

Now, we must add the constraints. Like design variables, these may be
added individually or by group. It is recommended that related
constraints are added by group where possible::

  >>> optProb.addConGroup('con',2, lower=None, upper=0.0)

This call adds two variables with name 'con'. There is no lower bound
for the variables and the upper bound is 0.0. 

We must also assign the the key value for the objective using the
``addObj()`` call::

  >>> optProb.addObj('obj')

The optimization problem can be printed to verify that it is setup correctly::

  >>> print optProb

To solve an optimization problem with ``pyOptSparse`` an optimizer
must be initialized. The initialization of one or more optimizers is
independent of the initialization of any number of optimization
problems. To initialize ``SLSQP``, which is an open-source, sequential
least squares programming algorithm that comes as part of the
pyOptSparse package, use::

  >>> opt = pyoptsparse.SLSQP()

This initializes an instance of ``SLSQP`` with the default options. The
setOption() method can be used to change any optimizer specific option,
for example the internal output flag of ``SLSQP``::

  >>> opt.setOption('IPRINT', -1)

Now TP37 can be solved using  ``SLSQP`` and for example, ``pyOptSparse``’s automatic
finite difference for the gradients::

  >>> sol = opt(optProb, sensType='FD')

We can print the solution objection to view the result of the optimization::

  >>> print sol

    TP037
  ================================================================================

          Objective Function: objfunc

      Solution: 
  --------------------------------------------------------------------------------
      Total Time:                    0.0256
         User Objective Time :       0.0003
         User Sensitivity Time :     0.0021
         Interface Time :            0.0226
         Opt Solver Time:            0.0007
      Calls to Objective Function :      23
      Calls to Sens Function :            9

      Objectives:
          Name        Value        Optimum
  	     f               0             0

  	  Variables (c - continuous, i - integer, d - discrete):
             Name      Type       Value       Lower Bound  Upper Bound
	    xvars_0     c	     24.000000       0.00e+00     4.20e+01 
	    xvars_1     c	     12.000000       0.00e+00     4.20e+01 
	    xvars_2     c	     12.000000       0.00e+00     4.20e+01 

   	  Constraints (i - inequality, e - equality):
          Name    Type                    Bounds
	      con   	  i        1.00e-20 <= 0.000000 <= 0.00e+00
	      con   	  i        1.00e-20 <= 0.000000 <= 0.00e+00

  --------------------------------------------------------------------------------
