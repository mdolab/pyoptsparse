.. _uno:

IPOPT
=====
`UNO <https://github.com/cvanaret/Uno>`_ (Unifying Nonlinear Optimization) is a C++ library with a Python interface.


Installation
------------
UNO and its Python interface must be installed via the unopy package  `pyuno <https://pypi.org/project/unopy/>`_.

Options
-------

UNO decomposes nonlinear optimization algoirthms into a set of "ingredents" that can be combined in various ways.
To mimick the behavior of well-known algorithms, it has several ``presets``.
Further customization is available via other options.

Please refer to the `UNO options documentation https://github.com/cvanaret/Uno/blob/main/docs/options.md` for a a complete list of options.

Please refer to the `IPOPT website <https://coin-or.github.io/Ipopt/OPTIONS.html>`_ for complete listing of options.
The following are the options which are set by default within pyOptSparse.
All other options take the default value with IPOPT and cyipopt unless specified by the user.

.. optionstable:: pyoptsparse.pyUNO.pyUNO.UNO
   :filename: IPOPT_options.yaml

Informs
-------
.. optionstable:: pyoptsparse.pyUNO.pyUNO.UNO
   :type: informs

API
---
.. currentmodule:: pyoptsparse.pyUNO.pyUNO

.. autoclass:: UNO
   :members: __call__
