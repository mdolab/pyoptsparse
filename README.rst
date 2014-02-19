pyOptSparse - PYthon OPTimization (Sparse) Framework
====================================================
Copyright (c) 2013-2014, Gaetan K. W. Kenway

pyOptSparse is fairly extensive modification to the original pyOpt
framework originally written by Dr. Ruben E. Perez and
Peter W. Jansen. It aims to be mostly backwards compatible with pyOpt;
existing optimization scripts can typicaly be converted to use
pyOptSparse with little to no modification. However, to use many of
the new features of pyOptSparse, some script modification will
typically be necesasry. 

The original goal of the pyOpt is to create an object-oriented
framework for formulating and solving nonlinear constrained
optimization problems.

From the pyOpt README some of the main pyOpt features are:

* Object-oriented development maintains independence between 
  the optimization problem formulation and its solution by 
  different optimizers
   
* Allows for easy integration of gradient-based, gradient-free, 
  and population-based optimization algorithms
    
* Interfaces both open source as well as industrial optimizers

* Ease the work required to do nested optimization and provides
  automated solution refinement

* On parallel systems it enables the use of optimizers when 
  running in a mpi parallel environment, allows for evaluation 
  of gradients in parallel, and can distribute function 
  evaluations for gradient-free optimizers

* Optimization solution histories can be stored during the 
  optimization process. A partial history can also be used 
  to warm-restart the optimization
    
After many years of experience of using pyOpt, the following list
summarizes the main modifications made to pyOpt to create pyOptSpase:

Changes to pyOptSparse
----------------------

* Elimination of n^2 scaling behaviour when adding large numbers of
  design variables (>10,000)

* Proper handling of all optimizers when run in parallel environment
  (Only a single optimization instance is run, not nProc instances)

* More flexible return specification of constraints
  
* More flexible return specification of constraint *gradients*
  
* Complete elimination of gradient indexing errors when using
  dictionary based returns

* Substantial improvement of optimization script robustness through
  indexing elimination
  
* Automatic assembly of sparse jacobians with both dense and sparse
  subblocks
  
* Automatic conversion of sparse jacobian to dense for optimizers that
  cannot handle sparse constraints

* User specification of design variable scaling
  
* User specificaiton of constraint scaling
    
* Time limit specification for many optimizers

* Default is now to use design variable groups
  
* Specification of linear constraints, dense or sparse (SNOPT only)

* Sparse non-linear jacobian (SNOPT only)
  
* New history file format. Uses naitve python shelve database
  format. Much easier to retrieve optimization data after
  optimization.

* Fixed hot start bug where first call to user functions is a
  gradient. It is now guaranteed, that the first call is to the
  function evlaution, not the gradient.

* Various bug fixes in SNOPT
  
* Equality constraints now work correctly across all optimizers
  
* Better optimizer class documentation with better formatting

* Python 3.x compatibility

pyOpt - Building and Installing
===============================
Copyright (c) 2013-2014, pyOpt Developers

Requirements
------------
pyOpt has the following dependencies:

* Python 2.5+ (but not Python 3.x)
* Numpy 1.0+
* Swig 1.3+
* c/FORTRAN compiler (compatible with f2py)
    
Further dependencies to take advantage 
of parallel computing capabilities are:

* mpi compiler
* mpi4py

Notes:

* In Windows MinGW is recommended if c/FORTRAN compilers are not available
* In Linux, the python header files (python-dev) are also required.
* Compatibility on Windows 64bit has not been tested 

Installation
------------
To install the pyOptSparse package in a folder on the Python search path 
(usually in a python site-packages or dist-packages folder) run:
    
>>> python setup.py install

This will typically require root access and thus the command actually needs to be:

>>> sudo python setup.py install

Alternatively, the recommended approcah to using pyOptSparse is to
install inplace. This does not require root access. 
    
>>> python setup.py inplace

To use pyOptSparse in this case, the user should add the path of the
diretory containing pyoptsparse to the user's PYTHONPATH enviromental
variable. For example, if the pyoptsparse directory is located at::

  /home/user/packages/pyoptsparse

The required line in the .bashrc file would be::

  export PYTHONPATH=$PYTHONPATH:/home/user/packages/

Notes:
    
    * You may want to uninstall any previous version of pyOpt before installing a new 
      version, as there may be conflicts.
    * Some optimizers are licensed and their sources are not included with this distribution. 
      To use them, please request their sources from the authors as indicated in the optimizer 
      LICENSE files, and place them in their respective source folders before installing the package.
    * In Windows, if MinGW is used make sure to install for it the C, C++, and Fortran compilers and run:
      
      >>> python setup.py install --compiler=mingw32
      
    * Installing to site-packages/ requires root privileges on Linux.
    * By default pyOpt will attempt to use compilers available on the system. To get a list of 
      available compilers and their corresponding flag on a specific system use:
      
      >>> python setup.py compilers

Tutorial
========

The following shows how to get started with pyOptSparse by solving
Schittkowski’s TP37 constrained problem. First, we show the complete
program listing and then go through each statement line by line::

  import pyoptsparse
  def objfunc(xx):
      x = xx['xvars']
      fobj = -x[0]*x[1]*x[2]
      conval = [0]*2
      conval[0] = x[0] + 2.*x[1] + 2.*x[2] - 72.0
      conval[1] = -x[0] - 2.*x[1] - 2.*x[2]
      fcon = {'con':conval}
      fail = False

      return fobj, fcon, fail

  optProb = pyoptsparse.Optimization('TP037', objfunc)
  optProb.addVarGroup('xvars',3, 'c',lower=[0,0,0], upper=[42,42,42], value=10)
  optProb.finalizeDesignVariables()
  optProb.addConGroup('con',2, lower=None, upper=0.0)
  print optProb
  opt = pyoptsparse.SLSQP()
  sol = opt(optProb, sens='FD')
  print sol

Start by importing the pyOptSparse package::

  >>> import pyoptsparse

Next we define the objective function that takes in the design
variable *dictionary* and returns the objective function value (a
scalar), a *dictionary* of constraints and a (boolean) flag indicating
if the objective function evaluation was successful. For the TP37, the
objective function is a simple analytic function::

  def objfunc(xdict):
      x = xdict['xvars']
      fobj = -x[0]*x[1]*x[2]
      conval = [0]*2
      conval0] = x[0] + 2.*x[1] + 2.*x[2] - 72.0
      conval[1] = -x[0] - 2.*x[1] - 2.*x[2]
      fcon = {'con':conval}
      fail = Flase

      return fobj, fcon, fail

Notes:

  1. The ``xdict`` variable is a dictionary whose keys are the names
     from each addVar() and addVarGroup() call. The line::

       x = xdict['xvars']

     retrieves an array of length 3 which are all the variables for
     this optimization. 

  2. The line::
    
       conval = [0]*2

     creates a list of length 2, which stores the numerical values of
     the two constraints. The constraint return, 'fcon' must be a
     dictionary whose keys much match the constraint names given in
     addCon() and addConGroup() calls.  This is done in the following call::
       
       fcon = {'con':conval}

Now the optimization problem can be initialized::

  >>> optProb = Optimization('TP037', objfunc)

This creates an instance of the optimization class with a name and a
reference to the objective function. To complete the setup of the
optimization problem, the design variables and constraints need to be defined. 

Design variables and constraints can be added either one-by-one or as
a group. Adding variables by group is generally recommended for
related variables::

  >>> optProb.addVarGroup('xvars', 3, 'c',lower=[0,0,0], upper=[42,42,42], value=10)

This calls adds a group of 3 variables with name 'xvars'. The variable
bounds (side constraints) are 0 for the lower bounds, 42 for the upper
bounds. The inital values for each variable is 10.0

After all variables have been added, the following function must be called::

  >>> optProb.finalizeDesignVariables()

After this call, no more variables may be added. Now, we must add the
constraints. Like design variables, these may be added individually or
by group. It is recommended that related constraints are added by
group where possible::

  >>> optProb.addConGroup('con',2, lower=None, upper=0.0)

This call adds two variables with name 'con'. There is no lower bound
for the variables and the upper bound is 0.0. 

The optimization problem can be printed to verify that it is setup correctly::

  >>> print optOpt

To solve an optimization problem with pyOptSparse an optimizer must be
initialized. The initialization of one or more optimizers is
independent of the initialization of any number of optimization
problems. To initialize SLSQP, which is an open-source, sequential
least squares programming algorithm that comes as part of the pyOptSparse
package, use::

  >>> opt = SLSQP()

This initializes an instance of SLSQP with the default options. The
setOption() method can be used to change any optimizer specific option,
for example the internal output flag of SLSQP::

  >>> slsqp.setOption('IPRINT', -1)

Now TP37 can be solved using SLSQP and for example, pyOptSparse’s automatic
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



Licensing
---------
Distributed using the GNU Lesser General Public License (LGPL); see 
the LICENSE file for details.

Please cite pyOpt and the authors of the respective optimization
algorithms in any publication for which you find it useful. 
(This is not a legal requirement, just a polite request.)

Contact and Feedback
--------------------
If you have questions, comments, problems, want to contribute to the
framework development, or want to report a bug, please contact the 
main developer:
    
* Gaetan K. W. Kenway (gaetank _ a t _ gmail.com)
