IPOPT
=====
IPOPT (Interior Point OPTimizer) is an open source interior point optimizer, designed for large-scale nonlinear optimization.
The source code can be found `here <https://www.coin-or.org/download/source/Ipopt/>`_.
The latest version we support is 3.13.2.

Installation
------------
IPOPT must be installed separately, then linked to pyOptSparse when building.
For the full installation instructions, please see `their documentation <https://coin-or.github.io/Ipopt/INSTALL.html>`_.
OpenMDAO also has a very helpful `script <https://github.com/OpenMDAO/build_pyoptsparse/>`_ which can be used to install IPOPT with other linear solvers.
Here we explain a basic setup using MUMPS as the linear solver, together with METIS adapted from the OpenMDAO script.

#. Download the tarball and extract it to ``$IPOPT_DIR`` which could be set to for example ``$HOME/packages/Ipopt``.

#. Install METIS, which can be used to improve the performance of the MUMPS linear solver.

   .. code-block:: bash

      # build METIS
      cd $IPOPT_DIR
      git clone https://github.com/coin-or-tools/ThirdParty-Metis.git
      cd ThirdParty-Metis
      ./get.Metis
      ./configure --prefix=$IPOPT_DIR
      make
      make install

#. Install MUMPS

   .. code-block:: bash

      # build MUMPS
      cd $IPOPT_DIR
      git clone https://github.com/coin-or-tools/ThirdParty-Mumps.git
      cd ThirdParty-Mumps
      ./get.Mumps
      ./configure --with-metis --with-metis-lflags="-L${IPOPT_DIR}/lib -lcoinmetis" \
           --with-metis-cflags="-I${IPOPT_DIR}/include -I${IPOPT_DIR}/include/coin-or -I${IPOPT_DIR}/include/coin-or/metis" \
           --prefix=$IPOPT_DIR CFLAGS="-I${IPOPT_DIR}/include -I${IPOPT_DIR}/include/coin-or -I${IPOPT_DIR}/include/coin-or/metis" \
           FCFLAGS="-I${IPOPT_DIR}/include -I${IPOPT_DIR}/include/coin-or -I${IPOPT_DIR}/include/coin-or/metis"
      make
      make install

#. Build IPOPT

   .. code-block:: bash

      # build IPOPT
      cd $IPOPT_DIR
      mkdir build
      cd build
      ../configure --prefix=${IPOPT_DIR} --disable-java --with-mumps --with-mumps-lflags="-L${IPOPT_DIR}/lib -lcoinmumps" \
           --with-mumps-cflags="-I${IPOPT_DIR}/include/coin-or/mumps"
      make
      make install     

#. You must add the IPOPT library path to the ``LD_LIBRARY_PATH`` variable for things to work right.
   This could be done for example by adding the following to your ``.bashrc``::

     export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$IPOPT_DIR/lib

   Furthermore, the environment variable ``$IPOPT_DIR`` must be set correctly in order to link to pyOptSparse.
   Alternatively, you can manually define the variables ``$IPOPT_LIB`` and ``$IPOPT_INC`` for the lib and include paths separately.


#. Now clean build pyOptSparse. Verify that IPOPT works by running the relevant tests.

.. note::

   To get IPOPT working with pyOptSparse when using another linear solver, several things must be changed.

   #. The ``setup.py`` file located in ``pyoptsparse/pyIPOPT`` must be updated accordingly.
      In particular, the ``libraries=`` line must be changed to reflect the alternate linear solver.
      For  example, for HSL you need to replace ``coinmumps`` and ``coinmetis`` with ``coinhsl``.
   #. The option ``linear_solver`` in the options dictionary must be changed.
      The default value can be changed in ``pyIPOPT.py`` so that this option does not need to be manually set in every run script.

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

