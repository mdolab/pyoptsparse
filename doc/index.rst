===========
pyOptSparse
===========
pyOptSparse is an object-oriented framework for formulating and solving nonlinear constrained optimization problems in an efficient, reusable, and portable manner.
Some key features of pyOptSparse include:

- Object-oriented development, which maintains independence between the optimization problem formulation and its solution by different optimizers
- The use of sparse matrices throughout the code, to more efficiently handle large-scale optimization problems
- Parallel model execution under MPI, both for expensive analyses that must be done in parallel, and for parallel function evaluation when using certain gradient-free optimizers
- The optimization histories can be stored during the optimization process, and a partial history can also be used to hot-restart the optimization
- A post-processing GUI utility called OptView to analyze optimization results


pyOptSparse is a fork of `pyOpt <http://www.pyopt.org/>`_.
However, it is not backwards compatible with pyOpt and thus optimization scripts will need to be modified to use pyOptSparse. 

Getting Started
===============
To get started, please see the :ref:`install` and the :ref:`quickstart`.

.. toctree::
   :maxdepth: 1
   :caption: Table of Contents
   :hidden:

   install
   quickstart
   guide
   advancedFeatures
   postprocessing
   changes
   contribute
   publishedWorks
   citation
   license

.. toctree::
   :maxdepth: 1
   :caption: API Reference
   :hidden:

   api/optimization
   api/optimizer
   api/gradient
   api/variable
   api/constraint
   api/objective
   api/solution
   api/history
   api/utils

.. toctree::
   :maxdepth: 1
   :caption: Optimizers
   :hidden:

   optimizers/SNOPT
   optimizers/IPOPT
   optimizers/SLSQP
   optimizers/NLPQLP
   optimizers/NSGA2
   optimizers/PSQP
   optimizers/ParOpt
   optimizers/CONMIN
   optimizers/ALPSO
