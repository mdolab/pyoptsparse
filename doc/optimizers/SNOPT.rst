.. _snopt:

SNOPT
=====
SNOPT is a sparse nonlinear optimizer that is particularly useful for
solving large-scale constrained problems with smooth objective
functions and constraints. The algorithm consists of a sequential
quadratic programming (SQP) algorithm that uses a smooth augmented
Lagrangian merit function, while making explicit provision for
infeasibility in the original problem and in the quadratic programming
subproblems. The Hessian of the Lagrangian is approximated using the
BFGS quasi-Newton update.

Installation
------------
SNOPT is available for purchase `here 
<http://www.sbsi-sol-optimize.com/asp/sol_snopt.htm>`_. Upon purchase, you should receive a zip file. Within the zip file, there is a folder called ``src``. To use SNOPT with pyoptsparse, paste all files from ``src`` except snopth.f into ``pyoptsparse/pySNOPT/source``.

From v2.0 onwards, only SNOPT v7.7.x is officially supported.
To use pyOptSparse with previous versions of SNOPT, please checkout release v1.2.
The current version of SNOPT being tested is v7.7.5, although v7.7.1 is also expected to work.


Options
-------
The default values come directly from SNOPT, which is why they are not displayed in the table below.
Please refer to the user manual for those options.
Note that there are two key differences:

- The SNOPT option ``Proximal iterations limit`` has its default value changed to 10000, in order to fully solve the proximal point problem to optimality
- The option ``Save major iteration variables`` is unique to the Python wrapper, and takes a list of values which can be saved at each iteration to the History file.
  Possible values are

  - ``step``
  - ``merit``
  - ``feasibility``
  - ``optimality``
  - ``penalty``
  - ``Hessian``
  - ``slack``
  - ``lambda``
  - ``condZHZ``

.. optimizertable:: pyoptsparse.pySNOPT.pySNOPT.SNOPT
   :type: options
   :nodefvalues:


Informs
-------
.. optimizertable:: pyoptsparse.pySNOPT.pySNOPT.SNOPT
   :type: informs

API
---
.. currentmodule:: pyoptsparse.pySNOPT.pySNOPT

.. autoclass:: SNOPT
   :members: __call__

