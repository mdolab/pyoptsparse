.. _paropt:

ParOpt
======
ParOpt is a nonlinear interior point optimizer that is designed for large parallel design optimization problems.
ParOpt is open source and can be downloaded at `https://github.com/smdogroup/paropt <https://github.com/smdogroup/paropt>`_.
Documentation and examples for ParOpt can be found at `https://smdogroup.github.io/paropt/ <https://smdogroup.github.io/paropt/>`_.
Unlike the other optimizers supported by pyOptSparse, the pyOptSparse interface to ParOpt is maintained as part of ParOpt itself.
We support ParOpt version 2.1.5 and later.

Installation
------------
Please follow the instructions `here <https://smdogroup.github.io/paropt/>`_ to install ParOpt as a separate Python package.
Make sure that the package is named ``paropt`` and the installation location can be found by Python, so that ``from paropt import ParOpt`` works within the pyOptSparse folder.

Options
-------
Please see the ParOpt documentation for all available options.

API
---
.. currentmodule:: pyoptsparse.pyParOpt.ParOpt

.. autoclass:: ParOpt
   :members: __call__
