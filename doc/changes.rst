Major changes compared to pyOpt
===============================

The following list summarizes some of the changes/improvements made to pyOpt to get pyOptSparse:

* Elimination of :math:`\mathcal{O} (n^2)` scaling behaviour when adding large numbers of
  design variables (>10,000)
* Proper handling of all optimizers when run in parallel environment (only a single optimization instance is run, not one instance per processor)
* More flexible return specification of constraints
* More flexible return specification of constraint *gradients*
* Complete elimination of gradient indexing errors using dictionary-based returns
* Substantial improvement of optimization script robustness through indexing elimination
* Automatic assembly of sparse Jacobians with both dense and sparse sub-blocks
* Automatic conversion of sparse Jacobian to dense for optimizers that cannot handle sparse constraints
* User specification of design variable scaling
* User specification of constraint scaling
* Design variable returns in dictionary format only
* Specification of linear constraints, dense or sparse 
* Sparse non-linear Jacobians 
* New history file format. Uses SQLite dictionaries.
* Fixed hot start bug where first call to user functions is a gradient. It is now guaranteed, that the first call is to the function evaluation, not the gradient.
* Various bug fixes in SNOPT
* Constraints can be in *any* order, independent of what the individual
  optimizer requires
* Python 3.x compatibility
