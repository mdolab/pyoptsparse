Advanced Features
=================
.. Gradient Evaluation with Complex Step
.. -------------------------------------

.. Parallel Execution
.. ------------------

.. Storing Optimization History
.. ----------------------------

Hot start
--------
Hot start refers to the way optimizations are initialized.
Suppose you run an optimization, and it was accidentally terminated prematurely, for example due to an iteration limit which was too low.
If you restarted the optimization from the last iteration, you will likely incur a performance penalty, and end up taking far more steps than it would've taken if the original optimization was allowed to progress.
This is because for many optimizers, the next iterate does not depend only on the current iterate, but effectively all the previous iterates as well.
In a sense, these optimizers are path-dependent, and when you restart the optimization you lose all the accumulated history.

Hot starting optimizations is a way to address this issue.
pyOptSparse has the ability to read in the previous optimization history file, and `replay` the entire history starting at the original design variables, feeding ``funcs`` to the optimizer each time.
For a deterministic optimizer, if all solver settings remain the same, then it should request the same sequence of iterates (up to machine precision), and if those have been previously evaluated, they are read from the history file and passed to the optimizer.
If at any point the requested design variables diverge from the history, then we will actually perform a function evaluation, and all subsequent points will be newly evaluated.
For this process to work, the following must be true:

-  The optimizer is deterministic.
   Note that a stochastic optimizer is still deterministic even if it uses random numbers, as long as the initial seed is fixed.
-  Optimizer settings that affect the path of the optimization must remain the same.
   It is perfectly fine to change for example settings related to print outs, but not those affecting the line search.

To use the hot start feature, simply call the optimizer with the option ``hotStart = <hot start file>``.
See the API documentation for each optimizer for more information.



.. Time limit
.. ----------

.. Clean Optimization Termination
.. ------------------------------
