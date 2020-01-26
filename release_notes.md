# Release Notes for pyOptSparse v1.0

October 22, 2019

v1.0 is the first major release of pyOptSparse.
As such, the release notes will only highlight a few recent developments, rather than serve as an exhaustive list.

## Bug Fixes: 
- various minor fixes to code testing

## New Features:
- ParOpt has been added as a new optimizer
- Optimal Lagrange multipliers are now saved in the solution object as `lambdastar`, for those optimizers that provide this variable

## Changes to code behavior:
- SNOPT will no longer perform the additional function call after finishing the optimization. This was a feature of SNOPT meant for the user to perform any finalization and clean up, and write out any necessary files prior to the end of the program. However, this option was never provided to pySNOPT users, and therefore has been disabled. This will reduce the number of function evaluations for every optimization by one, with no impact on the final result. **However, the user can no longer rely on using `db[db['last']]` to retrieve the optimal design, since it may not be the last function evaluation.** The optimal design vector is stored under the key `xs`.