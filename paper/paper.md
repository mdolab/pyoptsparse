---
title: 'pyOptSparse: A Python framework for large-scale constrained nonlinear optimization of sparse systems'
tags:
  - optimization
  - Python
authors:
  - name: Neil Wu
    orcid: 0000-0001-8856-9661
    affiliation: 1
  - name: Gaetan Kenway
    affiliation: 1
  - name: Charles A. Mader
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

pyOptSparse considers optimization problems of the form
\begin{align*}
\text{minimize}\quad & f(x)\\
\text{with respect to}\quad & x\\
\text{such that}\quad & l \le \begin{pmatrix}
x\\
Ax\\
g(x)\\
\end{pmatrix}
\le u\\
\end{align*}
where $x$ is the vector of design variables and $f(x)$ is a nonlinear objective function.
$A$ is the linear constraint Jacobian, and $g(x)$ is the set of nonlinear constraint functions.
pyOptSparse makes a distinction between linear and nonlinear constraints, since some optimizers have special treatments for linear constraints.
Sparse linear constraints can be directly supplied in a number of different formats, and pyOptSparse will automatically handle the assembly and conversion for each optimizer.
For nonlinear constraints, the sparsity pattern of the constraint Jacobian $\nabla g(x)$ can be specified as well.

# Features
## Support for multiple optimizers
pyOptSparse provides built-in support for a number of popular proprietary and open-source optimizers.
By isolating the optimization problem definition from the optimizer, the user is able to easily switch between different optimizers applied to the same optimization problem.
In fact, a single line of code change is required, with minimal effort from the user.

Although pyOptSparse is primarily focused on large-scale gradient-based optimization, it provides support for gradient-free optimizers as well.
In addition, discrete variables, multi-objective and population-based optimizers are all supported.
Due to the object-oriented programming approach, it is also straightforward to extend pyOptSparse to support any additional optimizers not currently available.
That way, new optimizers would automatically inherit all the capabilities within pyOptSparse.

## String-based indexing
Unlike many of the other optimization frameworks out there, pyOptSparse is designed for large-scape optimizations, particularly in engineering applications.
With up to $\mathcal{O}(1000)$ design variables and constraints, it is important to keep track of all values during computation.
pyOptSparse employs string-based indexing to accomplish this.
Instead of using a single flattened array, related design variables and constraints can be grouped together into arrays.
These are combined into a Python dictionary, with each group identified by a unique key.
Similarly, the constraint Jacobian is represented using a nested dictionary approach.
This has several advantages:
- Accessing design variable and constraint values without needing to consider their global indices, which is clumsy and prone to error.
  The global indices are also often optimizer-dependent, and this extra level of abstraction allows pyOptSparse to handle everything behind the scenes.
- The constraint Jacobian can be computed and provided at the sub-block level, with the whole Jacobian assembled by pyOptSparse.
  This closely mirrors the engineering workflow, where different sub-blocks of the Jacobian are often computed by different tools.
  The user only has to ensure that the indices within each sub-block is correct, and the rest is handled automatically.

## Support for sparse linear and nonlinear constraints
One major feature of pyOptSparse is the support for sparse constraints.
When defining constraints, it's possible to provide the sparsity pattern of the Jacobian.
This can be done at the global level, by specifying which constraints are independent of which design variables, thereby letting pyOptSparse know that the corresponding sub-blocks of the Jacobian are always zero.
For non-zero sub-blocks, it's also possible to supply the sparsity pattern of the Jacobian, again using local indexing, such that the actual derivative computation can use sparse matrices as well.

pyOptSparse also provides explicit support for linear constraints, since some optimizers can satisfy linear constraints exactly at each major iteration.
In these cases, only the Jacobian, upper and lower bounds of the constraint need to be supplied.
The values and gradients of these constraints do not need to be evaluated every iteration, since the optimizer will satisfy them internally.

## Automatic computation of derivatives
If analytic derivatives for the objective and constraint functions are not available, pyOptSparse can automatically compute them internally using finite differences or complex step.
For finite differences, the user has the option of either using forward or central differences, and with either an absolute or relative step size.
Computing derivatives using finite differences can be expensive, requiring $n+1$ evaluations for forward differences and $2n+1$ for centered differences.
They are also inaccurate due to subtractive cancellation error under finite precision arithmetic.
Complex step [@Martins2003] on the other hand, still requires the same number of function evaluations.
However, by perturbing each design variable along the imaginary axis, subtractive cancellation errors are avoided.
As a result, by taking very small steps, derivatives can be accurate down to machine precision.
Of course, attention must be paid to ensure that the objective and constraint functions can be evaluated correctly with complex design variables.

## Optimizer-independent problem scaling
pyOptSparse offers optimizer-independent scaling for individual design variables, objective and constraints.
By separating the optimization problem definition from the optimizer used, pyOptSparse is able to apply the scaling automatically and consistently to any chosen optimizer.
Since the optimization problem is always defined in the physical, user-defined space, bounds on the design variables and constraints do not need to be modified when applying a different scaling.
Furthermore, for gradient-based optimizers, all the derivatives are scaled automatically and consistently without any effort from the user.
The only thing to specify is a `scale` option when defining design variables, objective and constraints, and to change the scaling amounts to changing a single number.
This is particularly useful in engineering applications, where the physical quantities can sometimes cause undesirable problem scaling, and leading to poor optimization convergence.
pyOptSparse allows the user to adjust problem scaling for each design variable, constraint and objective separately, without needing to change the bound specification or derivative computation.

## Parallel execution
<!-- mention mpi4py? -->
pyOptSparse offers parallel execution in three ways.
Firstly and most commonly, it can perform parallel function evaluations when the functions themselves require multiple processors.
This is commonly the case when performing large-scale optimizations, where the objective and constraint functions are the result of a complex analysis such as computational fluid simulations.
In this scenario, pyOptSparse can be executed with multiple processors, and all processors will be used to perform function evaluation, but only the root processor will run the optimizer itself.
That way, we avoid the scenario where each processor is running an independent copy of the optimizer, and potentially causing inconsistencies or processor locking.

Secondly, it is possible to perform parallel gradient evaluation when automatic finite difference or complex step derivatives are used.
If the function evaluation only requires a single processor, it's possible to call pyOptSparse with multiple processors so that each point in the finite difference stencil is evaluated in parallel, reducing the wall time for derivative computation.

Lastly, certain population-based optimizers may support parallel function evaluation for each optimizer iteration.
In the case of genetic algorithm or particle swarm optimization, a fixed number of function evaluations are required at each optimizer iteration, and they can be done in parallel if multiple processors are available and the functions only require a single processor to execute.
However, the implementation of this mechanism is optimizer-dependent, and support varies between the different optimizers available.

## Optimization visualization and restart
pyOptSparse can store an optimization history file using its own format based on SQLite.
The optimization history can then be visualized using the OptView tool provided.
Alternatively, users can post-process results manually.
To facilitate this task, an API is provided for querying the history file and accessing the optimization history.

The history file also enables two types of optimization restarts.
A *cold start* can be used by simply setting the initial design variables to the final design variables of a previous optimization.
A *hot start*, on the other hand, will replay the previous optimization in its entirety.
For a deterministic optimizer, as long as the functions and gradients remain the same, it will generate the same sequence of iterates.
pyOptSparse will simply retrieve the previously-evaluated quantities and supply them to the optimizer, without calling the objective function.
This feature is particularly useful if the objective function is expensive to evaluate, and the optimization terminated due to problems such as reaching the maximum iteration limit.
In this case, the full state within the optimizer can be regenerated through the hot start process, so that the optimization can continue without performance penalties.

# Simple optimization script
To highlight some of the features discussed above, we present the pyOptSparse script to solve a toy problem involving two sets of design variables $x$ and $y$, subject to two nonlinear constraints and one linear constraint.
The optimization is as follows:
\begin{align*}
\text{minimize}\quad & x_0 + x_1^3 + y_0^2 + y_1^2 + y_2^2 + y_3^2\\
\text{with respect to}\quad & x_0, x_1, y_0, y_1, y_2, y_3\\
\text{such that}\quad & -10 \le x_0 x_1\\
& -10 \le 3x_0 - \sin(x_1) \le 10\\
& y_0 - 2y_1 = 5
\end{align*}

The sparsity structure of the constraint Jacobian is shown below:
```
                 x (2)   y (4)
               +---------------+
       con (2) |   X   |       |
               -----------------
lin_con(L) (1) |       |   X   |
               +---------------+
```
which allows us to only specify derivatives for the two nonzero sub-blocks.
For simplicity, we will supply the linear Jacobian explicitly, and use complex step to compute the derivative for the nonlinear constraints automatically.

First we define the imports and the objective function.
```python
import numpy as np
from pyoptsparse import Optimization, OPT

def objfunc(xdict):
    x = xdict["x"]
    y = xdict["y"]
    funcs = {}
    funcs["obj"] = x[0] + x[1] ** 3 + np.sum(np.power(y, 2))
    funcs["con"] = np.zeros(2, np.complex)
    funcs["con"][0] = x[0] * x[1]
    funcs["con"][1] = 3 * x[0] - np.sin(x[1])
    fail = False
    return funcs, fail
```
Note that only the nonlinear constraints need to be evaluated here.

Next, we set up the optimization problem, including design variables, objective and constraints.
```python
# Optimization Object
optProb = Optimization("Example Optimization", objfunc)

# Design Variables
nx = 2
lower = [-10, -10]
upper = [10, 100]
value = [-5, 6]
optProb.addVarGroup("x", nx, lower=lower, upper=upper, value=value)
ny = 4
optProb.addVarGroup("y", ny, lower=-10, upper=None, value=0)

# Nonlinear constraints
ncons = 2
lower = [-10, -10]
upper = [None, 10]
optProb.addConGroup("con", ncons, wrt="x", lower=lower, upper=upper)
# Linear constraint
jac = np.zeros((1, ny))
jac[0, 0] = 1
jac[0, 1] = -2
optProb.addConGroup("lin_con", 1, lower=5, upper=5, wrt=["y"], jac={"y": jac}, linear=True)
# Objective
optProb.addObj("obj")
```
By using the `wrt` argument when adding constraints, we tell pyOptSparse that only the specified sub-blocks of the Jacobian are nonzero.

The linear Jacobian for this problem is
\begin{align*}\begin{pmatrix}
1\\
-2\\
0\\
0
\end{pmatrix}\end{align*}
which we construct as `jac` and pass to pyOptSparse.
For large optimization problems, the Jacobian can be constructed using sparse matrices.

Finally, we set up the optimizer and solve the optimization problem.
```python
# Optimizer
optOptions = {}
opt = OPT("SLSQP", options=optOptions)

# Optimize
sol = opt(optProb, sens="CS")
print(sol)
```

# Statement of Need
pyOptSparse is a fork of pyOpt [@Perez2012], and as the name suggests, its primary motivation is to support the use of sparse linear and nonlinear constraints in the context of gradient-based optimization.
This sets itself apart from other optimization frameworks such as SciPy [@SciPy] and NLopt [@NLopt], which do not provide the same level of support for sparse constraints.
NLopt does not support sparse Jacobians, either for linear or nonlinear constraints.
While SciPy does support both, it does not allow for the Jacobians to be specified separately for each sub-block, since it treats the design vector as a single array.
These frameworks also do not offer convenience features such as user-supplied optimization problem scaling, optimization hot-start, or post-processing utilities.
Although pyOptSparse is a general optimization framework, it is tailored to gradient-based optimizations of large-scale problems with sparse constraints.

pyOptSparse has been used extensively in engineering applications, particularly in the field of multidisciplinary design optimization (MDO).
Researchers have used it to perform aerodynamic shape optimization of aircraft wings [@Secco2019] and wind turbines [@Madsen2019], and aero-structural optimization of an entire aircraft [@Brooks2018].
pyOptSparse is also supported by OpenMDAO [@Gray2019], a popular Python framework for multidisciplinary analysis and optimization.
Through OpenMDAO, pyOptSparse has been applied to problems such as low-fidelity aero-structural wing design [@Chauhan2020] and aeropropulsive optimization of a boundary-layer ingestion propulsor [@Gray2018].

# Acknowledgements
We acknowledge the efforts of the original pyOpt developers, particularly Dr. Ruben E. Perez and Dr. Peter W. Jansen, who laid the foundation of the code.
We also acknowledge the numerous pyOptSparse users who have contributed to the code over the years.

# References
