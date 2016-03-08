.. _ipopt:

IPOPT
-----

Install instructions for ``pyIPOPT``.

#. Download ``IPopt-3.11.7.tar.gz`` and put in the ``/pyoptsparse/pyoptsparse/pyIPOPT`` directory

#. Untar

#. Rename the directory from ``Ipopt-x.x.x`` to ``Ipopt``

#. Copy ``ma27ad.f`` from the ``pyOptSparse`` bitbucket page into the ``Ipopt/ThirdParty/HSLold/`` directory

#. Go to::

     Ipopt/ThirdParty/Blas/

   and run::

      sh ./get.Blas

   This will download a blas copy and ``Ipopt`` will use that.

#. Go to::

     Ipopt/ThirdParty/Lapack/

   and run::

      sh ./get.Lapack

#. Run in the root directory
   ::
     $ ./configure --disable-linear-solver-loader

#. Now make
   ::
     $ make install

#. You must add the ``lib`` directory ``Ipopt`` to your
   ``LD_LIBRARY_PATH`` variable for things to work right
   ::
     export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/user/hg/pyoptsparse/pyoptsparse/pyIPOPT/Ipopt/lib

#. Now the pyOptSparse builder (run from the root directory) should take care of the rest. 


.. currentmodule:: pyoptsparse.pyoptsparse.pyIPOPT.pyIPOPT

.. autoclass:: IPOPT
   :members: __call__

