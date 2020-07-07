---
title: 'pyOptSparse: A Python package for constrained nonlinear optimization of large sparse systems'
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
<!-- Add some general description of pyOptSparse, cite [@Perez2012a] -->
<!-- Cite some papers here that have used pyOptSparse -->

# The optimization problem
In pyOptSparse, we consider optimization problems of the form:
\begin{align}
\min\quad &f(x)\\
\text{with respect to}\quad &x\\
\text{such that}\quad g_{j,\text{L}} &\le g_j(x) \le g_{j,\text{U}}, \quad j = 1, ..., m\\
x_{i,\text{L}} &\le x_i \le x_{i,\text{U}}, \quad i = 1, ..., n
\end{align}
where $x$ is the vector of $n$ design variables, $f(x)$ is a nonlinear function, and $g(x)$ is a set of $m$ nonlinear functions.

Equality constraints are specified using the same upper and lower bounds for the constraint, i.e., $g_{j,\text{L}} = g_{j,\text{U}}$.

# Available optimizers
pyOptSparse provides Python wrappers for the following optimizers.
However, some optimizers have a restrictive license and is therefore not provided as source code with the pyOptSparse package.

## ALPSO
Augmented Lagrangian Particle Swarm Optimizer [@Jansen2011a] is a gradient-free optimizer that uses particle swarm optmization (PSO), together with an augmented Lagrange multiplier method to enforce constraints.
It also allows for parallel function evaluations.

## CONMIN
CONMIN [@Vanderplaats1973] uses the method of feasible directions to solve constrained problems, by selecting a feasible search direction and step size at each iteration which  improves the objective function.

## IPOPT
Interior Point OPTimizer [@Curtis:2010:A] is an interior-point method that implements a primal-dual filter line search algorithm for large-scale nonlinear optimization problems.

## NLPQLP
NLPQLP [@Schittkowski2006] is a variant of NLPQL, allowing for parallel function evaluations during the line search.
The algorithm itself uses a sequential quadratic programming (SQP) approach based on the Lagrangian and linearized constraints, together with a non-monotone line search using an augmented-Lagrangian merit function.

## NSGA2
Non-dominated Sorting Genetic Algorithm II [@Deb2002a] is an evolutionary algorithm that can be applied to both single and multi-objective problems, particularly those that are non-convex and non-smooth.
By emphasizing non-dominated populations, the algorithm can converge to the optimal Pareto front for multi-objective problems.

## ParOpt
ParOpt [@Chin2019] is an interior-point optimizer well-suited for topology and multimaterial optimization problems.
In these types of problems, a large number of sparse linear constraints are typically present, and ParOpt takes advantage of the sparsity structure of the constraint Jacobians to solve the resultant linear systems in parallel.

## PSQP
Preconditioned Sequential Quadratic Programming is an optimizer which uses SQP together with a BFGS variable metric update to solve constrained nonlinear optimization problems.

## SLSQP
SLSQP [@Kraft1988a] is a sequential least squares programming algorithm which uses the Han–Powell quasi–Newton method with a BFGS update and an $\ell_1$ test function in the line search.

## SNOPT
Sparse NOnlinear OPTimizer [@Gill2005a] is an SQP optimizer designed for large-scale nonlinear constrained problems.
It is specifically designed for problems with a smooth objective and a large number of sparse constraints.

# Features
In this section we highlight some salient features of pyOptSparse that distinguishes itself from other optimizer packages.

## Exploiting the problem sparsity

## Parallelism

## Optimization restart

## Post-processing results

# References
