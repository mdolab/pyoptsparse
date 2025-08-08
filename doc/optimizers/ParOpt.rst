.. _paropt:

ParOpt
======
ParOpt is an open source package that implements trust-region, interior points, and MMA optimization algorithms.
The ParOpt optimizers are themselves MPI parallel, which allows them to scale to large problems.
Unlike other optimizers supported by pyOptSparse, as of ParOpt version 2.1.5 and later, the pyOptSparse interface to ParOpt is a part of ParOpt itself.
Maintaining the wrapper, and control of which versions of pyOptSparse are compatible with which versions of ParOpt, is therefore the responsibility of the ParOpt developers.
ParOpt can be downloaded at `<https://github.com/smdogroup/paropt>`_.
Documentation and examples can be found at `<https://smdogroup.github.io/paropt/>`_.
The wrapper code in pyOptSparse is minimal, simply allowing the ParOpt wrapper to be used in the same way as other optimizers in pyOptSparse, through the ``OPT`` method.

The ParOpt wrapper takes a ``sparse`` argument, which controls whether ParOpt uses sparse or dense storage for the constraint Jacobian.
The default is ``True``, which uses sparse storage, but is incompatible with ParOpt's trust-region algorithm.
If you want to use the trust-region algorithm, you must set ``sparse=False``, e.g.:

.. code-block:: python

   from pyoptsparse import OPT
   opt = OPT("ParOpt", sparse=False)

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
