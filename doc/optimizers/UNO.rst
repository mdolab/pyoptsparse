.. _uno:

UNO
===

`UNO <https://github.com/cvanaret/Uno>`_ (Unifying Nonlinear Optimization) is a C++ library with a Python interface.

Installation
------------

UNO and its Python interface must be installed via the unopy package `unopy <https://pypi.org/project/unopy/>`_.

Options
-------

UNO decomposes nonlinear optimization algorithms into a set of "ingredients" that can be combined in various ways.
To mimic the behavior of well-known algorithms, it has several ``presets``.
Further customization is available via other options.

These options may not reflect all options available in ``unopy``.
For complete and up-to-date information about all options, please refer to the `UNO options documentation <https://github.com/cvanaret/Uno/blob/main/docs/options.md>`_.

.. optionstable:: pyoptsparse.pyUNO.pyUNO.UNO
   :filename: UNO_options.yaml

Informs
-------
.. optionstable:: pyoptsparse.pyUNO.pyUNO.UNO
   :type: informs

API
---
.. currentmodule:: pyoptsparse.pyUNO.pyUNO

.. autoclass:: UNO
   :members: __call__
