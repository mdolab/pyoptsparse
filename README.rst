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
