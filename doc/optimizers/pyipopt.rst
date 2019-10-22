.. _ipopt:

IPOPT
=====
IPOPT (Interior Point OPTimizer) is an open source interior point optimizer, designed for large-scale nonlinear optimization.
The source code can be found `here <https://github.com/coin-or/Ipopt>`_.
The latest version we support is 3.11.7.

.. warning:: Currently only Python 2 is supported for IPOPT. Python 3 support
     will be provided once the IPOPT support is upgraded to a more recent package.

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

API
---
.. currentmodule:: pyoptsparse.pyIPOPT.pyIPOPT

.. autoclass:: IPOPT
   :members: __call__

