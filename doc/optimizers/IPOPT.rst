.. _ipopt:

IPOPT
=====
IPOPT (Interior Point OPTimizer) is an open source interior point optimizer, designed for large-scale nonlinear optimization.
The source code can be found `here <https://www.coin-or.org/download/source/Ipopt/>`_.
The latest version we support is 3.14.17.

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
All other options take the default value with IPOPT unless specified by the user.

.. optionstable:: pyoptsparse.pyIPOPT.pyIPOPT.IPOPT
   :filename: IPOPT_options.yaml


Informs
-------
.. optionstable:: pyoptsparse.pyIPOPT.pyIPOPT.IPOPT
   :type: informs

API
---
.. currentmodule:: pyoptsparse.pyIPOPT.pyIPOPT

.. autoclass:: IPOPT
   :members: __call__

