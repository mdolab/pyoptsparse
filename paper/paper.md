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
    affiliation: 2
  - name: Joaquim R. R. A. Martins
    affiliation: 1
affiliations:
 - name: Department of Aerospace Engineering, University of Michigan
   index: 1
 - name: NREL
   index: 2
date: June 24, 2020
bibliography: paper.bib
---

# Summary
pyOptSparse is an optimization framework designed for constrained nonlinear optimization of large sparse problems, providing a unified interface for a variety of gradient-free and gradient-based optimizers.
By using an object-oriented approach, it maintains independence between the optimization problem formulation and the implementation of the specific optimizers.
The code is MPI-wrapped to enable execution of expensive parallel analyses and gradient evaluations, such as using computational fluid dynamics (CFD) simulations.
During optimization, a history file can be stored to record the optimization history, both for post-processing and for performing an optimization hot-start.
A graphical user interface application is also provided, which can interactively plot various quantities over an optimization.

pyOptSparse is a fork of pyOpt [@Perez2012], and as the name suggests, its primary motivation is to support the use of sparse linear and nonlinear constraints which sets itself apart from other frameworks such as SciPy [@SciPy] and NLopt [NLopt].
pyOptSparse considers optimization problems of the form
\begin{align}
\text{min}\quad & f(x)\\
\text{with respect to}\quad & x\\
\text{such that}\quad & l \le \begin{pmatrix}
x\\
Ax\\
g_j(x)\\
\end{pmatrix}
\le u\\
\end{align}
where $x$ is the vector of design variables and $f(x)$ is a nonlinear objective function.
$A$ is the linear constraint Jacobian, and $g(x)$ is the set of nonlinear constraint functions.
pyOptSparse makes a distinction between linear and nonlinear constraints, since some optimizers have special treatments for linear constraints.
Sparse linear constraints can be directly supplied in a number of different formats, and pyOptSparse will automatically handle the assembly and conversion for each optimizer.
For nonlinear constraints, the sparsity pattern of the constraint Jacobian $\nabla g(x)$ can be specified as well.

pyOptSparse has been used extensively in the field of multidisciplinary design optimization (MDO).
Researchers have used it to perform aerodynamic shape optimization of aircraft wings[@Secco2019] and wind turbines [@Madsen2019], and aerostructural optimization of an entire aircraft [@Brooks2018].
pyOptSparse is also supported by OpenMDAO [@Gray2019], a popular framework for multidisciplinary analysis and optimization.
Through OpenMDAO, pyOptSparse has been applied to problems such as low-fidelity aero-structural wing design [@Chauhan2020] and aeropropulsive optimization of a boundary-layer ingestion propulsor [@Gray2018].

<!-- # Available optimizers
pyOptSparse provides Python wrappers for the following optimizers.
However, some optimizers have restrictive licenses and are therefore not provided as source code with the pyOptSparse package.

**ALPSO** [@Jansen2011] is a gradient-free optimizer that uses particle swarm optimization (PSO), together with an augmented Lagrange multiplier method to enforce constraints.
It also allows for parallel function evaluations.

**CONMIN** [@Vanderplaats1973] uses the method of feasible directions to solve constrained problems, by selecting a feasible search direction and step size at each iteration which  improves the objective function.

**IPOPT** [@Curtis2010] is an interior-point method that implements a primal-dual filter line search algorithm for large-scale nonlinear optimization problems.

**NLPQLP** [@Schittkowski2006] is a variant of NLPQL, allowing for parallel function evaluations during the line search.
The algorithm itself uses a sequential quadratic programming (SQP) approach based on the Lagrangian and linearized constraints, together with a non-monotone line search using an augmented-Lagrangian merit function.

**NSGA2** [@Deb2002] is an evolutionary algorithm that can be applied to both single and multi-objective problems, particularly those that are non-convex and non-smooth.
By emphasizing non-dominated populations, the algorithm can converge to the optimal Pareto front for multi-objective problems.

**ParOpt** [@Chin2019] is an interior-point optimizer well-suited for topology and multimaterial optimization problems.
In these types of problems, a large number of sparse linear constraints are typically present, and ParOpt takes advantage of the sparsity structure of the constraint Jacobians to solve the resultant linear systems in parallel.

**PSQP** is an optimizer which uses SQP together with a BFGS variable metric update to solve constrained nonlinear optimization problems.

**SLSQP** [@Kraft1988] is a sequential least squares programming algorithm which uses the Han–Powell quasi–Newton method with a BFGS update and an $\ell_1$ test function in the line search.

**SNOPT** [@Gill2005] is an SQP optimizer designed for large-scale nonlinear constrained problems.
It is specifically designed for problems with a smooth objective and a large number of sparse constraints.

# Acknowledgements
pyOpt -->

# References
