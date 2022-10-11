## Solve test problem HS15 from the Hock & Schittkowski collection.
#
#  min   100 (x2 - x1^2)^2 + (1 - x1)^2
#  s.t.  x1 x2 >= 1
#        x1 + x2^2 >= 0
#        x1 <= 0.5
#
#  The standard start point (-2, 1) usually converges to the standard
#  minimum at (0.5, 2.0), with final objective = 306.5.
#  Sometimes the solver converges to another local minimum
#  at (-0.79212, -1.26243), with final objective = 360.4.
##

# Standard Python modules
import argparse

# External modules
import numpy as np

# First party modules
from pyoptsparse import OPT, Optimization

parser = argparse.ArgumentParser()
parser.add_argument("--opt", help="optimizer", type=str, default="SLSQP")
parser.add_argument("--storeHistory", help="option to store history", action="store_true")
args = parser.parse_args()
optOptions = {}


def objfunc(xdict):
    x = xdict["xvars"]
    funcs = {}
    funcs["obj"] = [100 * (x[1] - x[0] ** 2) ** 2 + (1 - x[0]) ** 2]
    conval = np.zeros(2, "D")
    conval[0] = x[0] * x[1]
    conval[1] = x[0] + x[1] ** 2
    funcs["con"] = conval
    fail = False

    return funcs, fail


def sens(xdict, funcs):
    x = xdict["xvars"]
    funcsSens = {
        "obj": {
            "xvars": [
                2 * 100 * (x[1] - x[0] ** 2) * (-2 * x[0]) - 2 * (1 - x[0]),
                2 * 100 * (x[1] - x[0] ** 2),
            ]
        },
        "con": {
            "xvars": [
                [x[1], x[0]],
                [1, 2 * x[1]],
            ]
        },
    }

    fail = False
    return funcsSens, fail


# Optimization Object
optProb = Optimization("HS15 Constraint Problem", objfunc)

# Design Variables
lower = [-5, -5]
upper = [0.5, 5]
value = [-2, 1]
optProb.addVarGroup("xvars", 2, lower=lower, upper=upper, value=value)

# Constraints
lower = [1, 0]
upper = [None, None]
optProb.addConGroup("con", 2, lower=lower, upper=upper)

# Objective
optProb.addObj("obj")

# Check optimization problem:
print(optProb)

# Optimizer
opt = OPT(args.opt, options=optOptions)

# Solution
if args.storeHistory:
    histFileName = f"{args.opt.lower()}_hs015_Hist.hst"
else:
    histFileName = None
# end
sol = opt(optProb, sens=sens, storeHistory=histFileName)

# Check Solution
print(sol)
