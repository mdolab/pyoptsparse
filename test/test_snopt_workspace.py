"""
This example is taken from the OpenOpt Examples website.

https://github.com/PythonCharmers/OOSuite/blob/master/FuncDesigner/FuncDesigner/examples/nlpSparse.py

The Jacobian is represented in dense format to test the workspace lengths in SNOPT.
"""

import unittest
import numpy as np
from numpy.testing import assert_allclose
from numpy import arange

from pyoptsparse import Optimization, OPT
from pyoptsparse.pyOpt_error import Error

N = 500


def objfunc(xdict):
    x = xdict["x"]
    y = xdict["y"]
    z = xdict["z"]
    funcs = {}
    funcs["obj"] = x ** 2 + 2 * np.sum(y ** 2) + 3 * np.sum(z)
    funcs["con1"] = x + 1e-3 * abs(x) ** 2.05
    funcs["con2"] = x ** 4 + np.sum(y) + np.sum(z ** 2)
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
            "x": 4 * x ** 3,
            "y": np.ones(N),
            "z": 2 * z,
        },
        "con3": {
            "x": 1.0,
            "z": np.ones(2 * N),
        },
    }

    return funcsSens, False


class TestSNOPTWorkspace(unittest.TestCase):
    def optimize(self, optName, optOptions={}, storeHistory=False):
        # Optimization Object
        optProb = Optimization("large and dense", objfunc)

        # Design Variables
        optProb.addVar("x", lower=-100, upper=150, value=0)
        optProb.addVarGroup("y", N, lower=-10 - arange(N), upper=arange(N), value=0)
        optProb.addVarGroup("z", 2 * N, upper=arange(2 * N), lower=-100 - arange(2 * N), value=0)

        # Constraints
        optProb.addCon("con1", upper=100, wrt=["x"])
        optProb.addCon("con2", upper=100)
        optProb.addCon("con3", lower=4, wrt=["x", "z"])
        optProb.addConGroup(
            "lincon",
            N,
            lower=2 - 3 * arange(N),
            linear=True,
            wrt=["x", "y"],
            jac={"x": np.ones((N, 1)), "y": np.eye(N)},
        )
        optProb.addObj("obj")

        # Optimizer
        try:
            opt = OPT(optName, options=optOptions)
        except Error:
            raise unittest.SkipTest("Optimizer not available:", optName)

        sol = opt(optProb, sens=sens)

        return sol

    def test_default(self):
        test_name = "SNOPT_workspace_default"
        optOptions = {
            "Print file": "{}.out".format(test_name),
            "Summary file": "{}_summary.out".format(test_name),
        }
        sol = self.optimize("snopt", optOptions=optOptions)

        # Check that overwriting the workspace lengths works
        tol = 1e-5
        assert_allclose(sol.objectives["obj"].value, 10.0, atol=tol, rtol=tol)
        assert_allclose(sol.variables["x"][0].value, 2.0, atol=tol, rtol=tol)

    def test_user(self):
        test_name = "SNOPT_workspace_user"
        optOptions = {
            "Print file": "{}.out".format(test_name),
            "Summary file": "{}_summary.out".format(test_name),
            "Total character workspace": 500,
        }
        sol = self.optimize("snopt", optOptions=optOptions)

        # Check that the workspace is too small without overwriting the lengths
        self.assertEqual(sol.optInform["value"], 84)


if __name__ == "__main__":
    unittest.main()
