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

7. Now the pyOptSparse builder should take care of the rest. 


.. currentmodule:: pyoptsparse.pyoptsparse.pyIPOPT

.. autoclass:: IPOPT
   :members: __call__

