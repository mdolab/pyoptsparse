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

* Elimination of n^2 scaling behaviour when adding large numbers of design variables (>10,000)

* Proper handling of all optimizers when run in parallel environment
  (Only a single optimization instance is run, not nProc instances)

* More flexible return specification of constraints
  
* More flexible return specification of constraint *gradients*
  
* Complete elimination of gradient indexing errors when using
  dictionary based returns

* Substantial improvement of optimization script robustness through
  indexing elimination
  
* Automatic assembly of sparse jacobians with both dense and sparse subblocks
  
* Automatic conversion of sparse jacobian to dense for optimizers that cannot
  handle sparse constraints

* User specification of design variable scaling.
  
* User specificaiton of constraint scaling
    
* Time limit specification for many optimizers

* Default is now to use design variable groups
  
* Specification of linear constraints (SNOPT only)

* Sparse non-linear jacobian (SNOPT only)
  
* New history file format. Uses naitve python shelve database
  format. Much easier to retrieve optimization data after
  optimization. 

* Various bug fixes in SNOPT
  
* Equality constraints now work correctly across all optimizers
  
* Better optimizer class documentation with better formatting

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
