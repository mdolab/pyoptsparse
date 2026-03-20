.. _uno:

Uno
===

`Uno <https://github.com/cvanaret/Uno>`_ (Unifying Nonlinear Optimization) is a C++ library for solving nonlinearly constrained optimization problems

Installation
------------

Uno and its Python interface must be installed via the unopy package `unopy <https://pypi.org/project/unopy/>`_.

Options
-------

Uno decomposes nonlinear optimization algorithms into a set of "ingredients" that can be combined in various ways.
To mimic the behavior of well-known algorithms, it has several ``presets``.
Further customization is available via other options.

For complete and up-to-date information about all options, please refer to the `Uno options documentation <https://github.com/cvanaret/Uno/blob/main/docs/options.md>`_.

.. optionstable:: pyoptsparse.pyUno.pyUno.Uno
   :filename: UNO_options.yaml

Informs
-------
.. optionstable:: pyoptsparse.pyUno.pyUno.Uno
   :type: informs

API
---
.. currentmodule:: pyoptsparse.pyUno.pyUno

.. autoclass:: Uno
   :members: __call__
