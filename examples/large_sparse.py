# This example is taken from the OpenOpt Examples website.
# http://trac.openopt.org/openopt/browser/PythonPackages/FuncDesigner/FuncDesigner/examples/nlpSparse.py

# Standard Python modules
import argparse

# External modules
import numpy as np
from scipy import sparse

# First party modules
from pyoptsparse import OPT, Optimization

N = 50000


def objfunc(xdict):
    x = xdict["x"]
    y = xdict["y"]
    z = xdict["z"]
    funcs = {}
    funcs["obj"] = x**2 + 2 * np.sum(y**2) + 3 * np.sum(z)
    funcs["con1"] = x + 1e-3 * abs(x) ** 2.05
    funcs["con2"] = x**4 + np.sum(y) + np.sum(z**2)
    funcs["con3"] = x + np.sum(z)

    return funcs, False


def sens(xdict, funcs):
    x = xdict["x"]
    y = xdict["y"]
    z = xdict["z"]
    funcsSens = {
        "obj": {
            "x": 2 * x,
            "y": 4 * y,
            "z": 3 * np.ones(2 * N),
        },
        "con1": {
            "x": 2.05 * x * (x * x) ** 0.025,
        },
        "con2": {
            "x": 4 * x**3,
            "y": np.ones(N),
            "z": 2 * z,
        },
        "con3": {
            "x": 1.0,
            "z": np.ones(2 * N),
        },
    }

    return funcsSens, False


def large_sparse(optimizer="SNOPT", optOptions=None):
    opt_options = {} if optOptions is None else optOptions

    # Optimization Object
    optProb = Optimization("large and sparse", objfunc)

    # Design Variables
    optProb.addVar("x", lower=-100, upper=150, value=0)
    optProb.addVarGroup("y", N, lower=-10 - np.arange(N), upper=np.arange(N), value=0)
    optProb.addVarGroup("z", 2 * N, upper=np.arange(2 * N), lower=-100 - np.arange(2 * N), value=0)
    # Constraints
    optProb.addCon("con1", upper=100, wrt=["x"])
    optProb.addCon("con2", upper=100)
    optProb.addCon("con3", lower=4, wrt=["x", "z"])
    optProb.addConGroup(
        "lincon",
        N,
        lower=2 - 3 * np.arange(N),
        linear=True,
        wrt=["x", "y"],
        jac={"x": np.ones((N, 1)), "y": sparse.spdiags(np.ones(N), 0, N, N)},
    )
    optProb.addObj("obj")
    # Optimizer
    opt = OPT(optimizer, options=opt_options)
    optProb.printSparsity()

    return opt, optProb


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--opt", help="optimizer", type=str, default="SNOPT")
    args = parser.parse_args()
    optOptions = {}
    opt, optProb = large_sparse(args.opt, optOptions)

    # Solution
    sol = opt(optProb, sens=sens)
    print(sol)
