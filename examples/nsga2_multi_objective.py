# First party modules
from pyoptsparse import NSGA2, Optimization


def objfunc(xdict):
    x = xdict["x"]
    y = xdict["y"]

    funcs = {}
    funcs["obj1"] = (x - 0.0) ** 2 + (y - 0.0) ** 2
    funcs["obj2"] = (x - 1.0) ** 2 + (y - 1.0) ** 2

    fail = False

    return funcs, fail


# Instantiate Optimization Problem
optProb = Optimization("Rosenbrock function", objfunc)
optProb.addVar("x", "c", value=0, lower=-600, upper=600)
optProb.addVar("y", "c", value=0, lower=-600, upper=600)

optProb.addObj("obj1")
optProb.addObj("obj2")

# 300 generations will find x=(0,0), 200 or less will find x=(1,1)
options = {"maxGen": 200}
opt = NSGA2(options=options)
sol = opt(optProb)

print(sol)
