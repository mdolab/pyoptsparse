.. _changes:

Changes to pyOptSparse
======================

The following list summarizes some of the changes/improvements made to
pyOpt to get ``pyOptSparse``:

* Elimination of n^2 scaling behaviour when adding large numbers of
  design variables (>10,000)

* Proper handling of all optimizers when run in parallel environment
  (Only a single optimization instance is run, not nProc instances)

* More flexible return specification of constraints
  
* More flexible return specification of constraint *gradients*
  
* Complete elimination of gradient indexing errors using dictionary
  based returns

* Substantial improvement of optimization script robustness through
  indexing elimination
  
* Automatic assembly of sparse jacobians with both dense and sparse
  subblocks
  
* Automatic conversion of sparse jacobian to dense for optimizers that
  cannot handle sparse constraints

* User specification of design variable scaling
  
* User specificaiton of constraint scaling
    
* Design variable returns in dictionary format only
  
* Specification of linear constraints, dense or sparse 

* Sparse non-linear jacobians 
  
* New history file format. Uses naitve python shelve database
  format. Much easier to retrieve optimization data after
  optimization.

* Fixed hot start bug where first call to user functions is a
  gradient. It is now guaranteed, that the first call is to the
  function evlaution, not the gradient.

* Various bug fixes in SNOPT
  
* Constraints can be *any* order, independent of what the individual
  optimizer requires
  
* Python 3.x compatibility
