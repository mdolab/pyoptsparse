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
The follow is a non-exhaustive list of works that have used pyOptSparse.

MACH-Aero
~~~~~~~~~
MACH-Aero is an open-source aerodynamic shape optimization framework developed by the MDO Lab at the University of Michigan.
It uses ADflow for Euler and RANS-based analysis and adjoint computation, and uses pyOptSparse for optimization.

.. bibliography::
   :list: bullet
   :filter: keywords % "MACH"

OpenMDAO
~~~~~~~~
OpenMDAO is a popular multidisciplinary design and optimization framework developed by NASA Glenn Research Center, and support pyOptSparse as one of two possible optimization drivers.

.. bibliography::
   :list: bullet
   :filter:  keywords % "OpenMDAO"


OpenAeroStruct
~~~~~~~~~~~~~~
OpenAeroStruct is a popular low-fidelity aerostructural optimization framework written in OpenMDAO, and as such can leverage the optimization capabilities of pyOptSparse.

.. bibliography::
   :list: bullet
   :filter:  keywords % "OpenAeroStruct"

SUAVE
~~~~~
SUAVE is a conceptual-level aircraft design framework developed at Stanford University that includes multiple analysis fidelities, and provides a pyOptSparse interface.

.. bibliography::
   :list: bullet
   :filter:  keywords % "SUAVE"

Other
~~~~~
.. bibliography::
   :list: bullet
   :filter:  keywords == ""