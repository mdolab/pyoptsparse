# Example of auto-refinement -- runs ALPSO followed by SNOPT.
# See TP109.py for more information.

# External modules
import numpy as np
from tp109 import objfunc  # Import objective from the other file

# First party modules
from pyoptsparse import ALPSO, SNOPT, Optimization

USE_LINEAR = True
# Optimization Object
optProb = Optimization("TP109 Constraint Problem", objfunc)

# Design Variables (Removed infinite bounds for ALPSO)
lower = [0.0, 0.0, -0.55, -0.55, 196, 196, 196, -400, -400]
upper = [2000, 2000, 0.55, 0.55, 252, 252, 252, 800, 800]
value = [0, 0, 0, 0, 0, 0, 0, 0, 0]
optProb.addVarGroup("xvars", 9, lower=lower, upper=upper, value=value)

# Constraints
lower = [0, 0, 0, 0, 0, 0, 0, 0]
upper = [None, None, 0, 0, 0, 0, 0, 0]
if not USE_LINEAR:
    lower.extend([0, 0])
    upper.extend([None, None])

optProb.addConGroup("con", len(lower), lower=lower, upper=upper)

# And the 2 linear constriants
if USE_LINEAR:
    jac = np.zeros((1, 9))
    jac[0, 3] = 1.0
    jac[0, 2] = -1.0
    optProb.addConGroup("lin_con", 1, lower=-0.55, upper=0.55, wrt=["xvars"], jac={"xvars": jac}, linear=True)

# Objective
optProb.addObj("obj")

# Check optimization problem:
print(optProb)
optProb.printSparsity()

# Global Optimizer: ALPSO
opt1 = ALPSO()

# Get first Solution
sol1 = opt1(optProb)
print(sol1)

# Now run the previous solution with SNOPT
opt2 = SNOPT()
sol2 = opt2(sol1)

# Check Solution
print(sol2)
