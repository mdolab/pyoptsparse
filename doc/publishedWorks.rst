pyOptSparse in published works
==============================
pyOptSparse has been used extensively in the field of engineering design optimization.
The following is a non-exhaustive list of works that have used pyOptSparse.
The citations are organized by the optimization framework for which pyOptSparse is used within.

MACH-Aero
---------
`MACH-Aero <https://github.com/mdolab/MACH-Aero>`__ is an open-source aerodynamic shape optimization framework developed by the MDO Lab at the University of Michigan.
Since the vast majority of publications from the MDO Lab uses the MACH framework (and as a result pyOptSparse), only a few select publications from the lab are listed below.
The rest are works in conjunction with collaborators.

.. bibliography::
   :list: bullet
   :filter: keywords % "MACH"

OpenMDAO
--------
`OpenMDAO <https://openmdao.org/>`_ is a popular multidisciplinary design and optimization framework developed by NASA Glenn Research Center, and support pyOptSparse as one of the possible optimization `drivers <https://openmdao.org/twodocs/versions/latest/features/building_blocks/drivers/index.html>`__.

.. bibliography::
   :list: bullet
   :filter:  keywords % "OpenMDAO"


OpenAeroStruct
~~~~~~~~~~~~~~
`OpenAeroStruct <https://github.com/mdolab/openaerostruct>`_ is a low-fidelity aerostructural optimization framework implemented in OpenMDAO, and as such can leverage the optimization capabilities of pyOptSparse.

.. bibliography::
   :list: bullet
   :filter:  keywords % "OpenAeroStruct"

SUAVE
-----
`SUAVE <https://suave.stanford.edu/>`_ is a conceptual-level aircraft design framework developed at Stanford University that includes multiple analysis fidelities, and provides a pyOptSparse interface.

.. bibliography::
   :list: bullet
   :filter:  keywords % "SUAVE"

Other
-----
The works which did not use any of the aforementioned frameworks are listed here.

.. bibliography::
   :list: bullet
   :filter:  keywords == ""
