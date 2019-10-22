.. _paropt:

ParOpt
======
ParOpt is a nonlinear interior point optimizer that is designed for large parallel design optimization problems with structured sparse constraints.
ParOpt is open source and can be downloaded at `https://github.com/gjkennedy/paropt <https://github.com/gjkennedy/paropt>`_.
Documentation and examples for ParOpt can be found at `https://gjkennedy.github.io/paropt/ <https://gjkennedy.github.io/paropt/>`_.
ParOpt does not provide version tagging, but the GitHub source code from September 2019 has been verified to work.

Installation
------------
Please follow the instructions `here <https://gjkennedy.github.io/paropt/>`_ to install ParOpt as a separate Python package.
Make sure that the package is named ``paropt`` and the installation location can be found by Python, so that ``from paropt import ParOpt`` works within the pyOptSparse folder.
This typically requires installing it in a location which is already present under ``$PYTHONPATH`` environment variable, or you can modify the ``.bashrc`` file and manually append the path.

API
---

.. currentmodule:: pyoptsparse.pyParOpt.ParOpt

.. autoclass:: ParOpt
   :members: __call__
