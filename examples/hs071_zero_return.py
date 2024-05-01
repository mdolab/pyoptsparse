# Standard Python modules
import argparse

# External modules
import numpy as np

# First party modules
from pyoptsparse import OPT, Optimization

parser = argparse.ArgumentParser()
parser.add_argument("--opt", help="optimizer", type=str, default="SLSQP")
args = parser.parse_args()
optOptions = {}


def objfunc(xdict):
    x0 = xdict["x0"][0]
    x1 = xdict["x1"][0]
    x2 = xdict["x2"][0]
    x3 = xdict["x3"][0]

    funcs = {}
    funcs["obj"] = x0 * x3 * (x0 + x1 + x2) + x2
    funcs["con1"] = [x0 * x1 * x2 * x3]
    funcs["con2"] = [x0 * x0 + x1 * x1 + x2 * x2 + x3 * x3]
    fail = False
    return funcs, fail


def sens(xdict, funcs):
    x0 = xdict["x0"][0]
    x1 = xdict["x1"][0]
    x2 = xdict["x2"][0]
    x3 = xdict["x3"][0]

    funcsSens = {}
    funcsSens["obj"] = {
        "x0": np.array([x0 * x3 + x3 * (x0 + x1 + x2)]),
        "x1": np.array([x0 * x3]),
        "x2": np.array([x0 * x3 + 1.0, 0]),
        "x3": np.array([x0 * (x0 + x1 + x2)]),
    }

    funcsSens["con1"] = {
        "x0": np.array([[x1 * x2 * x3]]),
        "x1": np.array([[x0 * x2 * x3]]),
        # 'x2': np.array([[x0*x1*x3, 0]]),
        "x3": np.array([[x0 * x1 * x2]]),
    }
    #    ^
    #    |
    # If we don't return any one of the constraint Jacobian blocks,
    # pyoptsparse will assume it to be zero.

    funcsSens["con2"] = {
        "x0": np.array([[2.0 * x0]]),
        "x1": np.array([[2.0 * x1]]),
        "x2": np.array([[2.0 * x2, 0]]),
        "x3": np.array([[2.0 * x3]]),
    }

    fail = False
    return funcsSens, fail


# Optimization Object
optProb = Optimization("HS071 Constraint Problem", objfunc)

# Design Variables
x0 = [1.0, 5.0, 5.0, 1.0]
optProb.addVarGroup("x0", 1, lower=1, upper=5, value=x0[0])
optProb.addVarGroup("x1", 1, lower=1, upper=5, value=x0[1])
optProb.addVarGroup("x2", 2, lower=1, upper=5, value=x0[2])
optProb.addVarGroup("x3", 1, lower=1, upper=5, value=x0[3])

# Constraints
optProb.addConGroup("con1", 1, lower=[25], upper=[1e19])
optProb.addConGroup("con2", 1, lower=[40], upper=[40])

# Objective
optProb.addObj("obj")

# Check optimization problem:
print(optProb)

# Optimizer
opt = OPT(args.opt, options=optOptions)

# Solution
sol = opt(optProb, sens=sens)

# Check Solution
print(sol)
