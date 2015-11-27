.. _guide:

Guide
-----

``pyOptSparse`` is designed to solve general, constrained nonlinear
optimization problems of the form::

  min f(x) w.r.t. x

  s.t. g_j_L <= g_j(x) <= g_j_U, j = 1, ..., m
  x_i_L <= x_i <= x_i_U, i = 1, ..., n

where:

``x`` is the vector of ``n`` design variables
``f(x)`` is a nonlinear function
``g(x)`` is a set of ``m`` nonlinear functions

Equality constraints are specified using the same upper and lower
bounds for the constraint. ie. g_j_L = g_j_U. The ordering of the
constraints is arbitrary; ``pyOptSparse`` reorders the problem
automatically depending on the requirements of each individual
optimizer.

The optimization class is created using the following call::

  >>> optProb = Optimization('name', objFun)

The general template of the objective function is as follows::

  def obj_fun(xdict):
    funcs = {}
    funcs['obj'] = function(x)
    funcs['con_name'] = function(x)
    fail = False # Or True if an analysis failed

    return funcs, fail

where:

 * ``funcs`` is the dictionary of constraints and objective value(s)

 * ``fail`` is a Boolean. False for successful evaluation and True for unsuccessful

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
keyword will result in the following::

  x_optimizer = x_user * scale

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

The simplest way to add a single constraint with no bounds (ie not a
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
bounds are given in their un-scaled form. Internally, ``pyOptSparse``
will use the scaling factor to produce the following constraint::

  con_optimizer = con_user * scale

In the example above, the constraint values are divided by 10000,
which results in a upper bound (that the optimizer sees) of 1.0. 

Constraints may also be flagged as liner using the ``linear=True``
keyword option. Some optimizers can perform special treatment on
linear constraint, often ensuring that they are always satisfied
exactly on every function call (SNOPT for example). Linear constraints
also require the use of the ``wrt`` and ``jac`` keyword
arguments. These are explained below. 

One of the major goals of ``pyOptSparse`` is to enable the use of
sparse constraint jacobians. (Hence the 'Sparse` in the name!).
Manually computing sparsity structure of the constraint Jacobian is
tedious at best and become even more complicated as optimization
scripts are modified by adding or deleting design variables and/or
constraints. ``pyOptSParse`` is designed to greatly facilitate the
assembly of sparse constraint jacobians, alleviating the user of thus
burden. The idea is that instead of the user computing a dense matrix
representing the constraint jacobian, a ``dictionary of keys``
approach is used which allows incrementally specifying parts of the
constraint jacobain. Consider the optimization problem given below::

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

The ``X``'s denote which parts of the jacobian have non-zero
values. ``pyOptSparse`` does not determine the sparsity structure of
the jacobian automatally, it must be specified by the user during
calls to ``addCon`` and ``addConGroup``.  By way of example, the code
that generates the  hypothetical optimization problem is as follows::

  optProb.addVarGroup('varA', 3)
  optProb.addVarGroup('varB', 1)
  optProb.addVarGroup('varC', 3)

  optProb.addConGroup('conA', 2, upper=0.0, wrt=['varB', 'varC'])
  optProb.addConGroup('conB', 2, upper=0.0, wrt=['varC', 'varA'])
  optProb.addConGroup('conC', 4, upper=0.0)
  optProb.addConGroup('conD', 3, upper=0.0, wrt=['varC'])

Note that the order of the ``wrt`` (which stands for with-respect-to)
is not significant. Furthermore, if the ``wrt`` argument is omitted
altogether, ``pyOptSparse`` assumes that the constraint is dense. 

Using the ``wrt`` keyword allows the user to determine the overall
sparsity structure of the constraint jacobian. However, we have
currently assumed that each of the blocks with an ``X`` in is a dense
sub-block. ``pyOptSparse`` allows each of the *sub-blocks* to itself
be sparse. ``pyOptSparse`` requires that this sparsity structure to be
specified when the constraint is added. This information is supplied
through the ``jac`` keyword argument. Lets say, that the (conD, varC)
block of the jacobian is actually a sparse and linear. By way of
example, the call instead may be as follows::

  jac = sparse.lil_matrix((3,3))
  jac[0,0] = 1.0
  jac[1,1] = 4.0
  jac[2,2] = 5.0

  optProb.addConGroup('conD', 3, upper=0.0, wrt=['varC'], linear=True, jac={'varC':jac})

We have created a linked list sparse matrix using
``scipy.sparse``. Any scipy sparse matrix format can be accepted. We
have then provided this constraint jacobian using the ``jac=`` keyword
argument. This argument is a dictionary, and the keys must match the
design variable sets given in the ``wrt`` to keyword. Essentially what
we have done is specified the which blocks of the constraint rows are
non-zero, and provided the sparsity structure of ones that are sparse. 

For linear constraints the values in ``jac`` are meaningful: They must
be the actual linear constraint jacobian values (which do not
change). For non-linear constraints, on the sparsity structure
(non-zero pattern) is significant. The values themselves will be
determined by a call the sens() function. 

Also note, that the ``wrt`` and ``jac`` keyword arguments are only
supported when user-supplied sensitivity is used. If one used the
automatic gradient in ``pyOptSparse`` the constraint jacobian will
necessarily be dense. 

Objectives
++++++++++

Each optimization will require at least one objective to be
added. This is accomplished using a the call::

  otpProb.addObj('obj')

What this does is tell ``pyOptSparse`` that the key ``obj`` in the
function returns will be taken as the objective. For optimizers that
can do multi-objective optimization, (NSGA2 for example) multiple
objectives can be added. Optimizers that can only handle one objective
enforce that only a single objective is added to the optimization description. 
