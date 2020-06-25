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

The optimization class is created using the following call::

  >>> optProb = Optimization('name', objFun)

The general template of the objective function is as follows:

.. code-block:: python

  def obj_fun(xdict):
    funcs = {}
    funcs['obj_name'] = function(xdict)
    funcs['con_name'] = function(xdict)
    fail = False # Or True if an analysis failed

    return funcs, fail

where:

 * ``funcs`` is the dictionary of constraints and objective value(s)

 * ``fail`` can be a Boolean or an int. False (or 0) for successful evaluation and True (or 1) for unsuccessful. Can also be 2 when using SNOPT and requesting a clean termination of the run.

If the Optimization problem is unconstrained, ``funcs`` will contain only the objective key(s).

Design Variables
++++++++++++++++

The simplest way to add a single continuous variable with no bounds
(side constraints) and initial value of 0.0 is::

   >>> optProb.addVar('var_name')

This will result in a scalar variable included in the ``x`` dictionary
call to ``obj_fun`` which can be accessed by doing::

  >>> x['var_name']

A more complex example will include lower bounds, upper bounds and a
non-zero initial value::

  >>> optProb.addVar('var_name',lower=-10, upper=5, value=-2)

The ``lower`` or ``upper`` keywords may be specified as ``None`` to
signify there is no bound on the variable.

Finally, an additional keyword argument ``scale`` can be specified
which will perform an internal design variable scaling. The ``scale``
keyword will result in the following:

.. math::

  x_\text{opt} = x_\text{user} \times \text{scale}

The purpose of the scale factor is ensure that design variables of
widely different magnitudes can be used in the same optimization. Is
it desirable to have the magnitude of all variables within an order of
magnitude or two of each other.

The ``addVarGroup`` call is similar to ``addVar`` except that it adds
a group of 1 or more variables. These variables are then returned as a
numpy array within the x-dictionary. For example, to add 10 variables
with no lower bound, and a scale factor of 0.1::

  >>> optProb.addVarGroup('con_group', 10, upper=2.5, scale=0.1)


Constraints
+++++++++++

The simplest way to add a single constraint with no bounds (i.e., not a
very useful constraint!) is::

  >>> optProb.addCon('not_a_real_constraint')

To include bounds on the constraints, use the ``lower`` and ``upper``
keyword arguments. If ``lower`` and ``upper`` are the same, it will be
treated as an equality constraint::

  >>> optProb.addCon('inequality_constraint', upper=10)
  >>> optProb.addCOn('equality_constraint', lower=5, upper=5)

Like design variables, it is often necessary to scale constraints such
that all constraint values are approximately the same order of
magnitude. This can be specified using the ``scale`` keyword::

  >>> optProb.addCon('scaled_constraint', upper=10000, scale=1.0/10000)

Even if the ``scale`` keyword is given, the ``lower`` and ``upper``
bounds are given in their un-scaled form. Internally, pyOptSparse
will use the scaling factor to produce the following constraint:

.. math::

  \text{con}_\text{opt} = \text{con}_\text{user} \times \text{scale}

In the example above, the constraint values are divided by 10000,
which results in a upper bound (that the optimizer sees) of 1.0.

Constraints may also be flagged as linear using the ``linear=True``
keyword option. Some optimizers can perform special treatment on
linear constraint, often ensuring that they are always satisfied
exactly on every function call (SNOPT for example). Linear constraints
also require the use of the ``wrt`` and ``jac`` keyword
arguments. These are explained below.

One of the major goals of pyOptSparse is to enable the use of
sparse constraint Jacobians. (Hence the 'Sparse' in the name!).
Manually computing sparsity structure of the constraint Jacobian is
tedious at best and become even more complicated as optimization
scripts are modified by adding or deleting design variables and/or
constraints. pyOptSparse is designed to greatly facilitate the
assembly of sparse constraint Jacobians, alleviating the user of this
burden. The idea is that instead of the user computing a dense matrix
representing the constraint Jacobian, a ``dictionary of keys``
approach is used which allows incrementally specifying parts of the
constraint Jacobian. Consider the optimization problem given below::

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

The ``X``'s denote which parts of the Jacobian have non-zero
values. pyOptSparse does not determine the sparsity structure of
the Jacobian automatically, it must be specified by the user during
calls to ``addCon`` and ``addConGroup``.  By way of example, the code
that generates the  hypothetical optimization problem is as follows:

.. code-block:: python

  optProb.addVarGroup('varA', 3)
  optProb.addVarGroup('varB', 1)
  optProb.addVarGroup('varC', 3)

  optProb.addConGroup('conA', 2, upper=0.0, wrt=['varB', 'varC'])
  optProb.addConGroup('conB', 2, upper=0.0, wrt=['varC', 'varA'])
  optProb.addConGroup('conC', 4, upper=0.0)
  optProb.addConGroup('conD', 3, upper=0.0, wrt=['varC'])

Note that the order of the ``wrt`` (which stands for with-respect-to)
is not significant. Furthermore, if the ``wrt`` argument is omitted
altogether, pyOptSparse assumes that the constraint is dense.

To examine the sparsity pattern, pyOptSparse can generate the ASCII table shown above.
To do so, use the following call after adding all the design variables, objectives and constraints::

  >>> optProb.printSparsity()

Using the ``wrt`` keyword allows the user to determine the overall
sparsity structure of the constraint Jacobian. However, we have
currently assumed that each of the blocks with an ``X`` in is a dense
sub-block. pyOptSparse allows each of the *sub-blocks* to itself
be sparse. pyOptSparse requires that this sparsity structure to be
specified when the constraint is added. This information is supplied
through the ``jac`` keyword argument. Lets say, that the (conD, varC)
block of the Jacobian is actually a sparse and linear. By way of
example, the call instead may be as follows:

.. code-block:: python

  jac = sparse.lil_matrix((3,3))
  jac[0,0] = 1.0
  jac[1,1] = 4.0
  jac[2,2] = 5.0

  optProb.addConGroup('conD', 3, upper=0.0, wrt=['varC'], linear=True, jac={'varC':jac})

We have created a linked list sparse matrix using
``scipy.sparse``. Any scipy sparse matrix format can be accepted. We
have then provided this constraint Jacobian using the ``jac`` keyword
argument. This argument is a dictionary, and the keys must match the
design variable sets given in the ``wrt`` to keyword. Essentially what
we have done is specified the which blocks of the constraint rows are
non-zero, and provided the sparsity structure of ones that are sparse.

For linear constraints the values in ``jac`` are meaningful: they must
be the actual linear constraint Jacobian values (which do not
change). For non-linear constraints, only the sparsity structure
(i.e. which entries are nonzero) is significant. The values themselves will be
determined by a call to the sens() function.

Also note, that the ``wrt`` and ``jac`` keyword arguments are only
supported when user-supplied sensitivity is used. If automatic gradients
from pyOptSparse are used, the constraint Jacobian will
necessarily be dense.

.. note::
    Currently, only the optimizers SNOPT and IPOPT support sparse Jacobians.

Objectives
++++++++++

Each optimization will require at least one objective to be
added. This is accomplished using a the call::

  optProb.addObj('obj_name')

What this does is tell pyOptSparse that the key ``obj_name`` in the
function returns will be taken as the objective. For optimizers that
can do multi-objective optimization (e.g. NSGA2), multiple
objectives can be added. Optimizers that can only handle one objective
enforce that only a single objective is added to the optimization description.
