.. _ipopt:

IPOPT
=====
`IPOPT <https://coin-or.github.io/Ipopt/>`_ (Interior Point OPTimizer) is an open source interior point optimizer, designed for large-scale nonlinear optimization.
The source code can be found `here <https://github.com/coin-or/Ipopt>`_.

Installation
------------
IPOPT and its Python interface `cyipopt <https://github.com/mechmotum/cyipopt>` must be installed separately.
Follow the instructions `here <https://cyipopt.readthedocs.io/en/stable/install.html>`_.
OpenMDAO also has a very helpful `script <https://github.com/OpenMDAO/build_pyoptsparse/>`_ which can be used to install IPOPT with other linear solvers,
but it does not install ``cyipopt`` for you.

Options
-------
Please refer to the `IPOPT website <https://coin-or.github.io/Ipopt/OPTIONS.html>`__ for complete listing of options.
The following are the options which are set by default within pyOptSparse.
All other options take the default value with IPOPT and cyipopt unless specified by the user.

.. optionstable:: pyoptsparse.pyIPOPT.pyIPOPT.IPOPT
   :filename: IPOPT_options.yaml

.. note::
  There are several options that significantly affect the IPOPT performance, and IPOPT's default values are not always the best depending on the problem characteristics.
  Here are several noteworthy options based on our experience, and users are encouraged to explore non-default values if needed.
  pyOptSparse does not override the default value for these options.
  * `nlp_scaling_method`:
  by default, IPOPT internally applies the scaling based on the gradient at the initial point.
  If the problem is already well scaled at pyOptSparse level, you may set this to `none` to disable IPOPT scaling.
  * `hessian_approximation`:
  since pyOptSparse does not support the Hessian callback yet, cyipopt automatically sets this to `limited-memory`.
  pyOptSparse users do not need to set this option manually.
  * `limited_memory_max_history`:
  this determines the number of most recent iterations that are used for the Hessian approximation.
  IPOPT's default is 6, but it is often better to set this to a larger value so the Hessian approximation can utilize more information.
  * `mu_init`:
  this is the initial value for the barrier parameter, and IPOPT's default is 0.1.
  This parameter has a significant impact on the search path, especially when you have a lot of constraints.
  If the initial point is good (i.e., feasible or close to feasible), setting a smaller value (e.g., 1e-5) often accelerates convergence.
  * `mu_strategy`
  This controls the strategy for updating the barrier parameter, and the IPOPT's default is `monotone`.
  The other option is `adaptive`, which may accelerate the convergence, but `monotone` tends to be more robust.


Informs
-------
.. optionstable:: pyoptsparse.pyIPOPT.pyIPOPT.IPOPT
   :type: informs

API
---
.. currentmodule:: pyoptsparse.pyIPOPT.pyIPOPT

.. autoclass:: IPOPT
   :members: __call__
