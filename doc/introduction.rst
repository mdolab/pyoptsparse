.. _introduction:

Introduction
============

``pyOptSparse`` is fairly extensive modification to the original pyOpt
framework originally written by Dr. Ruben E. Perez and
Peter W. Jansen. It is not backwards compatible with pyOpt and thus
optimization scripts will need to be modified to use ``pyOptSparse``. 

The original goal of pyOpt was to create an object-oriented
framework for formulating and solving nonlinear constrained
optimization problems. This goal remains the same with ``pyOptSparse``.

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
