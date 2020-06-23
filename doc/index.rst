===========
pyOptSparse
===========
pyOptSparse is an object-oriented framework for formulating and solving nonlinear constrained optimization problems in an efficient, reusable, and portable manner.
It is a fork of pyOpt that uses sparse matrices throughout the code to more efficiently handle large-scale optimization problems.

.. toctree::
   :maxdepth: 1

   introduction
   changes
   install
   tutorial
   guide
   postprocessing
   contribute

.. toctree::
   :maxdepth: 1
   :caption: API Reference
   :hidden:

   api/optimization
   api/optimizer
   api/constraint
   api/variable
   api/gradient
   api/history

.. toctree::
   :maxdepth: 1
   :caption: Optimizers
   :hidden:

   optimizers/pysnopt
   optimizers/pyipopt
   optimizers/pyslsqp
   optimizers/pynlpqlp
   optimizers/pynsga2
   optimizers/pypsqp
   optimizers/pyparopt
   optimizers/pyconmin
   optimizers/pyalpso