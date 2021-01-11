.. _citation:

Citation
========
If you use pyOptSparse, please cite the following paper:

    N. Wu, G. Kenway, C. A. Mader, J. Jasa, and J. R. R. A. Martins. pyOptSparse: A Python framework for large-scale constrained nonlinear optimization of sparse systems. Journal of Open Source Software, 5(54), 2564, October 2020. https://doi.org/10.21105/joss.02564

The paper is available online from the Journal of Open Source Software `here <https://joss.theoj.org/papers/10.21105/joss.02564>`__.
To cite this paper, you can use the following BibTeX entry:

.. code-block:: bibtex

    @article{Wu2020,
        doi = {10.21105/joss.02564},
        year = {2020},
        publisher = {The Open Journal},
        volume = {5},
        number = {54},
        pages = {2564},
        author = {Neil Wu and Gaetan Kenway and Charles A. Mader and John Jasa and Joaquim R. R. A. Martins},
        title = {pyOptSparse: A Python framework for large-scale constrained nonlinear optimization of sparse systems},
        journal = {Journal of Open Source Software}
    }


Applications
------------
pyOptSparse has been used extensively in the field of engineering design optimization.
The follow is a non-exhaustive list of works that have used pyOptSparse.
The citations are organized by the optimization framework for which pyOptSparse is used within.

MACH-Aero
~~~~~~~~~
`MACH-Aero <https://github.com/mdolab/MACH-Aero>`_ is an open-source aerodynamic shape optimization framework developed by the MDO Lab at the University of Michigan.
Since the vast majority of publications from the MDO Lab uses the MACH framework (and as a result pyOptSparse), only a few select publications from the lab are listed below.
The rest are works in conjunction with collaborators.

.. bibliography::
   :list: bullet
   :filter: keywords % "MACH"

OpenMDAO
~~~~~~~~
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
~~~~~
`SUAVE <https://suave.stanford.edu/>`_ is a conceptual-level aircraft design framework developed at Stanford University that includes multiple analysis fidelities, and provides a pyOptSparse interface.

.. bibliography::
   :list: bullet
   :filter:  keywords % "SUAVE"

Other
~~~~~
The works which did not use any of the aforementioned frameworks are listed here.

.. bibliography::
   :list: bullet
   :filter:  keywords == ""