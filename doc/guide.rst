Guide
-----

pyOptSparse is designed to solve general, constrained nonlinear
optimization problems of the form:

.. math::
  \min\quad &f(x)\\
  \text{with respect to}\quad &x\\
  \text{such that}\quad g_{j,\text{L}} &\le g_j(x) \le g_{j,\text{U}}, \quad j = 1, ..., m\\
  x_{i,\text{L}} &\le x_i \le x_{i,\text{U}}, \quad i = 1, ..., n

where:
:math:`x` is the vector of :math:`n` design variables,
:math:`f(x)` is a nonlinear function,
and :math:`g(x)` is a set of :math:`m` nonlinear functions.

Equality constraints are specified using the same upper and lower bounds for the constraint. i.e., :math:`g_{j,\text{L}} = g_{j,\text{U}}`.
The ordering of the constraints is arbitrary; pyOptSparse reorders the problem automatically depending on the requirements of each individual optimizer.

The optimization class is created using the following call:

.. code-block:: python

  optProb = Optimization("name", objconFun)

The general template of the objective and constraint function is as follows:

.. code-block:: python

  def objconFun(xdict):
      funcs = {}
      funcs["obj_name"] = function(xdict)
      funcs["con_name"] = function(xdict)
      fail = False  # Or True if an analysis failed

      return funcs, fail

where:

 * ``funcs`` is the dictionary of constraints and objective value(s)

 * ``fail`` can be a Boolean or an int. False (or 0) for successful evaluation and True (or 1) for unsuccessful. Can also be 2 when using SNOPT and requesting a clean termination of the run.

If the Optimization problem is unconstrained, ``funcs`` will contain only the objective key(s).

Design Variables
++++++++++++++++
The simplest way to add a single continuous variable with no bounds (side constraints) and initial value of 0.0 is
to simply call :meth:`addVar <pyoptsparse.pyOpt_optimization.Optimization.addVar>`:

.. code-block:: python

   optProb.addVar("var_name")

This will result in a scalar variable included in the ``x`` dictionary call to ``obj_fun`` which can be accessed by doing

.. code-block:: python

  x["var_name"]

A more complex example will include lower bounds, upper bounds and a non-zero initial value:

.. code-block:: python

  optProb.addVar("var_name", lower=-10, upper=5, value=-2)

The ``lower`` or ``upper`` keywords may be specified as ``None`` to signify there is no bound on the variable.

Finally, an additional keyword argument ``scale`` can be specified which will perform an internal design variable scaling.
The ``scale`` keyword will result in the following:

.. math::

  x_\text{opt} = x_\text{user} \times \text{scale}

The purpose of the scale factor is ensure that design variables of widely different magnitudes can be used in the same optimization.
It is desirable to have the magnitude of all variables within an order of magnitude or two of each other.

The :meth:`addVarGroup <pyoptsparse.pyOpt_optimization.Optimization.addVarGroup>` call is similar to 
:meth:`addVar <pyoptsparse.pyOpt_optimization.Optimization.addVar>` except that it adds a group of 1 or more variables.
These variables are then returned as a numpy array within the x-dictionary.
For example, to add 10 variables with no lower bound, and a scale factor of 0.1:

.. code-block:: python

  optProb.addVarGroup("con_group", 10, upper=2.5, scale=0.1)


Constraints
+++++++++++

The simplest way to add a single constraint with no bounds (i.e., not a very useful constraint!) is
to use the function :meth:`addCon <pyoptsparse.pyOpt_optimization.Optimization.addCon>`:

.. code-block:: python

  optProb.addCon("not_a_real_constraint")

To include bounds on the constraints, use the ``lower`` and ``upper`` keyword arguments.
If ``lower`` and ``upper`` are the same, it will be treated as an equality constraint:

.. code-block:: python

  optProb.addCon("inequality_constraint", upper=10)
  optProb.addCon("equality_constraint", lower=5, upper=5)

Like design variables, it is often necessary to scale constraints such that all constraint values are approximately the same order of magnitude.
This can be specified using the ``scale`` keyword:

.. code-block:: python

  optProb.addCon("scaled_constraint", upper=10000, scale=1.0 / 10000)

Even if the ``scale`` keyword is given, the ``lower`` and ``upper`` bounds are given in their un-scaled form.
Internally, pyOptSparse will use the scaling factor to produce the following constraint:

.. math::

  \text{con}_\text{opt} = \text{con}_\text{user} \times \text{scale}

In the example above, the constraint values are divided by 10000, which results in a upper bound (that the optimizer sees) of 1.0.

Constraints may also be flagged as linear using the ``linear=True`` keyword option.
Some optimizers can perform special treatment on linear constraint, often ensuring that they are always satisfied
exactly on every function call (SNOPT for example).
Linear constraints also require the use of the ``wrt`` and ``jac`` keyword arguments.
These are explained below.

One of the major goals of pyOptSparse is to enable the use of sparse constraint Jacobians, hence the `Sparse` in the name!
Manually computing sparsity structure of the constraint Jacobian is tedious at best and become even more complicated
as optimization scripts are modified by adding or deleting design variables and/or constraints.
pyOptSparse is designed to greatly facilitate the assembly of sparse constraint Jacobians, alleviating the user of this burden.
The idea is that instead of the user computing a dense matrix representing the constraint Jacobian,
a "dictionary of keys" approach is used which allows incrementally specifying parts of the constraint Jacobian.
Consider the optimization problem given below::

              varA (3)   varB (1)   varC (3)
            +--------------------------------+
   conA (2) |          |     X    |     X    |
            ----------------------------------
   conB (2) |     X    |          |     X    |
            ----------------------------------
   conC (4) |     X    |     X    |     X    |
            ----------------------------------
   conD (3) |          |          |     X    |
            +--------------------------------+

The ``X``'s denote which parts of the Jacobian have non-zero values.
pyOptSparse does not determine the sparsity structure of the Jacobian automatically,
it must be specified by the user during calls to :meth:`addCon <pyoptsparse.pyOpt_optimization.Optimization.addCon>` and :meth:`addConGroup <pyoptsparse.pyOpt_optimization.Optimization.addConGroup>`.
By way of example, the code that generates the  hypothetical optimization problem is as follows:

.. code-block:: python

  optProb.addVarGroup("varA", 3)
  optProb.addVarGroup("varB", 1)
  optProb.addVarGroup("varC", 3)

  optProb.addConGroup("conA", 2, upper=0.0, wrt=["varB", "varC"])
  optProb.addConGroup("conB", 2, upper=0.0, wrt=["varC", "varA"])
  optProb.addConGroup("conC", 4, upper=0.0)
  optProb.addConGroup("conD", 3, upper=0.0, wrt=["varC"])

Note that the order of the ``wrt`` (which stands for with-respect-to) is not significant.
Furthermore, if the ``wrt`` argument is omitted altogether, pyOptSparse assumes that the constraint is dense.

To examine the sparsity pattern, pyOptSparse can generate the ASCII table shown above.
To do so, use the following call after adding all the design variables, objectives and constraints:

.. code-block:: python

  optProb.printSparsity()

Using the ``wrt`` keyword allows the user to determine the overall sparsity structure of the constraint Jacobian.
However, we have currently assumed that each of the blocks with an ``X`` in is a dense sub-block.
pyOptSparse allows each of the *sub-blocks* to itself be sparse.
pyOptSparse requires this sparsity structure to be specified when the constraint is added.
This information is supplied through the ``jac`` keyword argument.
Lets say, that the ``(conD, varC)`` block of the Jacobian is actually a sparse and linear.
By way of example, the call instead may be as follows:

.. code-block:: python

  jac = sparse.lil_matrix((3, 3))
  jac[0, 0] = 1.0
  jac[1, 1] = 4.0
  jac[2, 2] = 5.0

  optProb.addConGroup("conD", 3, upper=0.0, wrt=["varC"], linear=True, jac={"varC": jac})

We have created a linked list sparse matrix using ``scipy.sparse``.
Any SciPy sparse matrix format can be accepted.
We have then provided this constraint Jacobian using the ``jac`` keyword argument.
This argument is a dictionary, and the keys must match the design variable sets given in the ``wrt`` to keyword.
Essentially what we have done is specified the which blocks of the constraint rows are non-zero,
and provided the sparsity structure of ones that are sparse.

Note that the ``wrt`` and ``jac`` keyword arguments are only supported when user-supplied sensitivity is used.
If automatic gradients from pyOptSparse are used, the constraint Jacobian will necessarily be dense.

.. note::
    Currently, only the optimizers SNOPT and IPOPT support sparse Jacobians.

Linear Constraints
~~~~~~~~~~~~~~~~~~
Linear constraints in pyOptSparse are defined exclusively by ``jac``, ``lower``, and ``upper`` entries of the ``addConGroup`` method.
For linear constraint :math:`g_L \leq Ax + b \leq g_U`, the constraint definition would look like:

.. code-block:: python
  
  optProb.addConGroup("con", num_cons, linear=True, wrt=["xvars"], jac={"xvars": A}, lower=gL - b, upper=gU - b)

Users should not provide the linear constraint values (i.e., :math:`g = Ax + b`) in a user-defined objective/constraint function.
pyOptSparse will raise an error if you do so.

For linear constraints, the values in ``jac`` are meaningful:
they must be the actual linear constraint Jacobian values (which do not change).
For non-linear constraints, only the sparsity structure (i.e. which entries are nonzero) is significant.
The values themselves will be determined by a call to the ``sens()`` function.


Objectives
++++++++++

Each optimization will require at least one objective to be added.
This is accomplished using a the call to :meth:`addObj <pyoptsparse.pyOpt_optimization.Optimization.addObj>`:

.. code-block:: python

  optProb.addObj("obj_name")

What this does is tell pyOptSparse that the key ``obj_name`` in the function returns will be taken as the objective.
For optimizers that can do multi-objective optimization (e.g. NSGA2), multiple objectives can be added.
Optimizers that can only handle one objective enforce that only a single objective is added to the optimization description.

Specifying Derivatives
++++++++++++++++++++++
Approximating Derivatives
~~~~~~~~~~~~~~~~~~~~~~~~~
pyOptSparse can automatically compute derivatives of the objective and constraint functions using finite differences or the complex-step method.
This is done by simply passing a string to the ``sens=`` argument when calling an optimizer.
See the possible values :ref:`here <gradient>`.
In the simplest case, using ``sens="FD"`` will be enough to run an optimization using forward differences with a default step size.

Analytic Derivatives
~~~~~~~~~~~~~~~~~~~~
If analytic derivatives are available, users can compute them within a user-defined function.
This function accepts as inputs a dictionary containing design variable values as well as another dictionary containing objective and constraint values.
It returns a nested dictionary containing the gradients of the objective and constraint values with respect to those design variables at the current design point.
Specifically, the first-layer keys should be associated with objective and constraint names while the second-layer keys correspond to design variables.
The dictionary values are the computed analytic derivatives, either in the form of lists or NumPy arrays with the expected shape.
Since pyOptSparse uses string indexing, users need to make sure the keys in the returned dictionary are consistent with the names of design variables, constraints and objectives which were first added to the optimization problem.

.. tip::
  #. Only the non-zero sub-blocks of the Jacobian need to be defined in the dictionary, and pyOptSparse will assume the rest to be zero.
  #. Derivatives of the linear constraints do not need to be given here, since they are constant and should have already been specified via the ``jac=`` keyword argument when adding the constraint.

For example, if the optimization problem has one objective ``obj``, two constraints ``con``, and three design variables ``xvars``, the returned sensitivity dictionary (with placeholder values) should have the following structure:

.. code-block:: python

  {"obj": {"xvars": [1, 2, 3]}, "con": {"xvars": [[4, 5, 6], [7, 8, 9]]}}

Once this function is constructed, users can pass its function handle to the optimizer when it's called via:

.. code-block:: python

  sol = opt(optProb, sens=sens, ...)


Optimizer Instantiation
+++++++++++++++++++++++
There are two ways to instantiate the optimizer object.
The first, and most explicit approach is to directly import the optimizer class, for example via

.. code-block:: python

  from pyoptsparse import SLSQP

  opt = SLSQP(...)

However, in order to easily switch between different optimizers without having to import each class, a convenience function called
:meth:`OPT <pyoptsparse.pyOpt_optimizer.OPT>` is provided.
It accepts a string argument in addition to the usual options, and instantiates the optimizer object based on the string:

.. code-block:: python

  from pyoptsparse import OPT

  opt = OPT("SLSQP", ...)

Note that the name of the optimizer is case-insensitive, so ``slsqp`` can also be used.
This makes it easy to for example choose the optimizer from the command-line, or more generally select the optimizer using strings without preemptively importing all classes.
