===========
pyOptSparse
===========
pyOptSparse is an object-oriented framework for formulating and solving nonlinear constrained optimization problems in an efficient, reusable, and portable manner.
Some key features of pyOptSparse include:

- Object-oriented development maintains independence between the optimization problem formulation and its solution by different optimizers
- Use of sparse matrices throughout the code to more efficiently handle large-scale optimization problems
- Enable the use of optimizers when running in a MPI parallel environment, allows for evaluation of gradients in parallel, and can distribute function evaluations for gradient-free optimizers
- Optimization solution histories can be stored during the optimization process. A partial history can also be used to warm-restart the optimization for any supported optimizer
- A post-processing GUI utility called OptView to analyze optimization results


pyOptSparse is a fork of `pyOpt <http://www.pyopt.org/>`_.
However, it is not backwards compatible with pyOpt and thus optimization scripts will need to be modified to use pyOptSparse. 

.. toctree::
   :maxdepth: 1
   :caption: Table of Contents
   :hidden:

   install
   quickstart
   guide
   postprocessing
   changes
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