---
title: 'pyOptSparse: A Python package for large-scale constrained nonlinear optimization of sparse systems'
tags:
  - optimization
  - Python
authors:
  - name: Neil Wu
    orcid: 0000-0001-8856-9661
    affiliation: 1
  - name: Gaetan Kenway
    affiliation: 1
  - name: John Jasa
    affiliation: 1
  - name: Joaquim R. R. A. Martins
    affiliation: 1
affiliations:
 - name: Department of Aerospace Engineering, University of Michigan
   index: 1
date: July 8, 2020
bibliography: paper.bib
---

# Summary
pyOptSparse is an optimization framework designed for constrained nonlinear optimization of large sparse problems, providing a unified interface for a variety of gradient-free and gradient-based optimizers.
By using an object-oriented approach, it maintains independence between the optimization problem formulation and the implementation of the specific optimizers.
The code is MPI-wrapped to enable execution of expensive parallel analyses and gradient evaluations, such as when using computational fluid dynamics (CFD) simulations which can require hundreds of processors.
During the optimization, a history file can be stored to record the optimization history, which can be used both for post-processing and for performing an optimization hot-start.
A graphical user interface application is also provided to interactively plot various quantities over an optimization.

pyOptSparse is a fork of pyOpt [@Perez2012], and as the name suggests, its primary motivation is to support the use of sparse linear and nonlinear constraints.
This sets itself apart from other optimization frameworks such as SciPy [@SciPy] and NLopt [@NLopt], and is more tailored to gradient-based optimizations of large-scale problems with sparse constraints.
pyOptSparse considers optimization problems of the form
\begin{align}
\text{minimize}\quad & f(x)\\
\text{with respect to}\quad & x\\
\text{such that}\quad & l \le \begin{pmatrix}
x\\
Ax\\
g(x)\\
\end{pmatrix}
\le u\\
\end{align}
where $x$ is the vector of design variables and $f(x)$ is a nonlinear objective function.
$A$ is the linear constraint Jacobian, and $g(x)$ is the set of nonlinear constraint functions.
pyOptSparse makes a distinction between linear and nonlinear constraints, since some optimizers have special treatments for linear constraints.
Sparse linear constraints can be directly supplied in a number of different formats, and pyOptSparse will automatically handle the assembly and conversion for each optimizer.
For nonlinear constraints, the sparsity pattern of the constraint Jacobian $\nabla g(x)$ can be specified as well.

pyOptSparse has been used extensively in the field of multidisciplinary design optimization (MDO).
Researchers have used it to perform aerodynamic shape optimization of aircraft wings [@Secco2019] and wind turbines [@Madsen2019], and aero-structural optimization of an entire aircraft [@Brooks2018].
pyOptSparse is also supported by OpenMDAO [@Gray2019], a popular Python framework for multidisciplinary analysis and optimization.
Through OpenMDAO, pyOptSparse has been applied to problems such as low-fidelity aero-structural wing design [@Chauhan2020] and aeropropulsive optimization of a boundary-layer ingestion propulsor [@Gray2018].

<!-- # Acknowledgements
pyOpt -->

# References
