.. _ipopt:

IPOPT
-----

Install instructions for ``pyOPT``.

1. Download ``IPopt-3.11.7.tar.gz`` and put in the ``/pyoptsparse/pyoptsparse/pyIPOPT`` directory

2. Untar

3. Rename the directory from ``Ipopt-x.x.x`` to ``Ipopt``

4. Copy ``ma27ad.f`` into the ``Ipopt/ThirdParty/HSLold/`` directory

5. Run in the root directory
   ::
     $ ./configure --disable-linear-solver-loader

6. Now make
   ::
     $ make install

7. You must add the ``lib`` directory ``Ipopt`` to your
   ``LD_LIBRARY_PATH`` variable for things to work right
   ::
     export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/user/hg/pyoptsparse/pyoptsparse/pyIPOPT/Ipopt/lib

8. Now the pyOptSparse builder should take care of the rest. 


.. currentmodule:: pyoptsparse.pyoptsparse.pyIPOPT

.. autoclass:: IPOPT
   :members: __call__

