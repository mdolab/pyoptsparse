#!/usr/bin/env python
from pyoptsparse import Optimization, OPT
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--sens", help="sensitivity mode", type=str, default="FD", choices=["FD", "CS", "CD", "user"])
parser.add_argument("--constrained", help="add a constraint", action="store_true")
parser.add_argument("--testHotstart", help="test hotstart", action="store_true")
parser.add_argument("--opt", help="optimizer", type=str, default="SLSQP")
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
opt = OPT(args.opt, options=optOptions)

if args.testHotstart:
    histName = "opt_hist.hst"
    # we map the max iteration option to the specific name for each optimizer
    if args.opt.lower() == "slsqp":
        OPTION = "MAXIT"
    elif args.opt.lower() == "snopt":
        OPTION = "Major iterations limit"
    elif args.opt.lower() == "ipopt":
        OPTION = "max_iter"
    elif args.opt.lower() == "ipopt":
        OPTION = "max_iter"
    elif args.opt.lower() == "nlpqlp":
        OPTION = "maxIt"
    elif args.opt.lower() == "psqp":
        OPTION = "MIT"
    elif args.opt.lower() == "paropt":
        OPTION = "tr_max_iterations"
    elif args.opt.lower() == "conmin":
        OPTION = "ITMAX"

    # First call just does 10 iterations
    opt.setOption(OPTION, 10)
    sol1 = opt(optProb, sens=sens, storeHistory=histName)

    # Now we are allowed to do 50
    opt.setOption(OPTION, 50)
    sol2 = opt(optProb, sens=sens, hotStart=histName, storeHistory=histName)
    print(sol2)
else:
    # Just do a normal run
    sol = opt(optProb, sens=sens)
    print(sol)
