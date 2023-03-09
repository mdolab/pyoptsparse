Advanced Features
=================
.. Gradient Evaluation with Complex Step
.. -------------------------------------

.. Parallel Execution
.. ------------------

MPI handling
------------
pyOptSparse can optionally run in parallel if a suitable ``mpi4py`` installation exists.
This will be automatically detected and imported at run-time.

If you only want to run in parallel, you can force pyOptSparse to do so by setting the environment variable
``PYOPTSPARSE_REQUIRE_MPI`` to any one of these values: ``['always', '1', 'true', 'yes']``
If a suitable ``mpi4py`` is not available, an exception will be raised and the run terminated.

If you explicitly do not wish to use ``mpi4py``, set the environment variable ``PYOPTSPARSE_REQUIRE_MPI`` to anything other than those values.
This can come in handy, for example, if your ``MPI`` installation is not functioning properly, but you still need to run serial code.

Storing Optimization History
----------------------------
pyOptSparse includes a :ref:`history` class that stores all the relevant optimization information an SQL database.
This database is updated at every optimization iteration, and can be accessed via both the API described in the linked section, and via :ref:`optview`.
By default, the history file is NOT written.
To turn the history recording on, use the ``storeHistory`` attribute when invoking the optimization run, e.g.:

.. code-block:: python

  sol = opt(optProb, sens=sens, storeHistory="<your-history-file-name>.hst", ...)


Hot start
---------
Hot start refers to the way optimizations are initialized.
Such start-up procedure is named in contrast to the regular optimization initialization, or "cold start", in which case pyOptSparse simply initializes an optimization job from scratch, using initial design variables set by the user.
There are several situations in which such cold starts can be avoided by leveraging on the information from a previous optimization, with the aim to reduce the overall computational time.
Suppose you run an optimization that was accidentally terminated prematurely by, for example, an excessively-low iteration limit.
If you restarted the optimization using the DVs from the last iteration (but losing all the accumulated history), you will start from a better initial point but will likely incur in a performance penalty.
The overall optimization, now split between two jobs, will end up taking far more steps than it would've taken if the original optimization was allowed to progress.
This is due to the fact that some of the optimizers wrapped in pyOptSparse need a few initial iterations to build a (local) approximation of the design space.
Because of this lack of information, there is a start-up cost to a cold started optimization, where the optimizer has to spend time rebuilding the information.
In a sense, for many optimizers the next iterate does not depend only on the current iterate, but effectively all the previous iterates as well!

Hot starting optimizations is a way to address this issue.
pyOptSparse has the ability to read in the previous optimization history file, and `replay` the entire history starting at the original design variables, feeding ``funcs`` to the optimizer each time.
For a deterministic optimizer, if all solver settings remain the same, then it should request the same sequence of iterates (up to machine precision), and if those have been previously evaluated, they are read from the history file and passed to the optimizer.
If at any point the requested design variables diverge from the history, then we will actually perform a function evaluation, and all subsequent points will be newly evaluated.
For this process to work, the following must be true:

-  The optimizer is deterministic.
   Note that a stochastic optimizer can still be made deterministic even if it uses random numbers, as long as the random seed is fixed.
   If the optimizer performs a random gradient check (e.g. SNOPT), it's best to disable these just in case.
-  Optimizer settings that affect the path of the optimization must remain the same.
   It is perfectly fine to change for example settings related to print outs, but not those affecting the line search.

To use the hot start feature, simply call the optimizer with the option ``hotStart = <hot start file>``.
See the API documentation for each optimizer for more information.
Because the hot start process will store all the previous "restarted" iterations in the new history file, it's possible to restart as many times as you like, each time using the previous history file.



Time limit (for SNOPT only)
---------------------------

.. note::

   Since SNOPT 7.7, the user should rely on the ``Time Limit`` SNOPT option instead of the ``timeLimit`` argument of the ``Optimizer`` instance to set the maximum optimization time.


The :ref:`optimizer` class in pyOptSparse has an attribute used to set the maximum allowable wall time for optimizations using SNOPT.
The code will exit gracefully when such time limit is reached.
This feature is particularly useful when running a time-constrained job, as in the case of most HPC systems.
To enable this feature, use the ``timeLimit`` option when invoking the optimizer, as shown below:

.. code-block:: python

  sol = opt(optProb, sens=sens, timeLimit=24 * 3600, ...)

Note that the attribute takes the maximum wall time *in seconds* as an integer number.

.. note::

   pyOptSparse will verify that the computational time is not exceeded before proceeding to the next iteration.
   It will NOT interrupt an ongoing function or sensitivity evaluation.
   If your function evaluations are expensive, you should be more conservative when setting the ``timeLimit`` option for it to be effective.


.. Clean Optimization Termination
.. ------------------------------
