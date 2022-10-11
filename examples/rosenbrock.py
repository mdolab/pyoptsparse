# Standard Python modules
import argparse

# First party modules
from pyoptsparse import SLSQP, Optimization

parser = argparse.ArgumentParser()
parser.add_argument("--sens", help="sensitivity mode", type=str, default="FD", choices=["FD", "CS", "CD", "user"])
parser.add_argument("--constrained", help="add a constraint", action="store_true")
parser.add_argument("--testHotstart", help="test hotstart", action="store_true")
args = parser.parse_args()
optOptions = {}


def objfunc(xdict):
    x = xdict["xvars"]  # Extract array
    funcs = {}
    funcs["obj"] = 100 * (x[1] - x[0] ** 2) ** 2 + (1 - x[0]) ** 2
    funcs["con"] = 0.1 - (x[0] - 1) ** 3 - (x[1] - 1)
    fail = False

    return funcs, fail


def sensfunc(xdict, funcsDict):
    x = xdict["xvars"]  # Extract array
    funcsSens = {}
    funcsSens["obj"] = {
        "xvars": [2 * 100 * (x[1] - x[0] ** 2) * (-2 * x[0]) - 2 * (1 - x[0]), 2 * 100 * (x[1] - x[0] ** 2)]
    }
    funcsSens["con"] = {"xvars": [-3 * (x[0] - 1) ** 2, -1]}
    fail = False
    return funcsSens, fail


if args.sens == "user":
    sens = sensfunc
else:
    sens = args.sens

# Instantiate Optimization Problem
optProb = Optimization("Rosenbrock function", objfunc)
optProb.addVarGroup("xvars", 2, "c", value=[3, -3], lower=-5.12, upper=5.12, scale=[1.0, 1.0])
if args.constrained:
    optProb.addCon("con", upper=0, scale=1.0)
optProb.addObj("obj")

# Create optimizer
opt = SLSQP(options=optOptions)

if args.testHotstart:
    histName = "opt_hist.hst"
    # First call just does 10 iterations
    opt.setOption("MAXIT", 10)
    sol1 = opt(optProb, sens=sens, storeHistory=histName)

    # Now we are allowed to do 50
    opt.setOption("MAXIT", 50)
    sol2 = opt(optProb, sens=sens, hotStart=histName, storeHistory=histName)
    print(sol2.fStar)
else:
    # Just do a normal run
    sol = opt(optProb, sens=sens)
    print(sol.fStar)
