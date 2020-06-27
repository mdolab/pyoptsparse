.. _ipopt:

IPOPT
=====
IPOPT (Interior Point OPTimizer) is an open source interior point optimizer, designed for large-scale nonlinear optimization.
The source code can be found `here <https://www.coin-or.org/download/source/Ipopt/>`_.
The latest version we support is 3.11.7.

Installation
------------

Install instructions for ``pyIPOPT``.

#. Download ``IPopt-3.11.7.tar.gz`` and put in the ``/pyoptsparse/pyoptsparse/pyIPOPT`` directory

#. Untar

#. Rename the directory from ``Ipopt-x.x.x`` to ``Ipopt``

#. Obtain the MA27 linear solver from `HSL <http://www.hsl.rl.ac.uk/download/MA27/1.0.0/a/>`_. Rename the source file ``ma27ad.f`` and put it in the ``Ipopt/ThirdParty/HSLold/`` directory

#. Go to::

     Ipopt/ThirdParty/Blas/

   and run::

     sh ./get.Blas

   This will download a blas copy and ``Ipopt`` will use that.

#. Go to::

     Ipopt/ThirdParty/Lapack/

   and run::

     sh ./get.Lapack

#. Run in the root directory::

     $ ./configure --disable-linear-solver-loader

#. Now make::

     $ make install

#. You must add the ``lib`` directory ``Ipopt`` to your
   ``LD_LIBRARY_PATH`` variable for things to work right::

     export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/user/hg/pyoptsparse/pyoptsparse/pyIPOPT/Ipopt/lib

#. Now the pyOptSparse builder (run from the root directory) should take care of the rest. 

Options
-------
.. optimizertable:: pyoptsparse.pyIPOPT.pyIPOPT.IPOPT
   :type: options

Informs
-------
.. optimizertable:: pyoptsparse.pyIPOPT.pyIPOPT.IPOPT
   :type: informs

API
---
.. currentmodule:: pyoptsparse.pyIPOPT.pyIPOPT

.. autoclass:: IPOPT
   :members: __call__

