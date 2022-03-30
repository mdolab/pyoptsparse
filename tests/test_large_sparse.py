"""
This example is taken from the OpenOpt Examples website.

https://github.com/PythonCharmers/OOSuite/blob/master/FuncDesigner/FuncDesigner/examples/nlpSparse.py

Only testing with SNOPT, which supports sparse formats.
"""

# Standard Python modules
import unittest

# External modules
import numpy as np
from parameterized import parameterized

# First party modules
from pyoptsparse import Optimization

# Local modules
from testing_utils import OptTest


class TestLarge(OptTest):
    name = "large_sparse"
    DVs = {"x", "y", "z"}
    objs = {"obj"}
    cons = {"con1", "con2", "con3"}
    xStar = {"x": 2}
    fStar = 10.0

    def objfunc(self, xdict):
        x = xdict["x"]
        y = xdict["y"]
        z = xdict["z"]
        funcs = {}
        funcs["obj"] = x**2 + 2 * np.sum(y**2) + 3 * np.sum(z)
        funcs["con1"] = x + 1e-3 * abs(x) ** 2.05
        funcs["con2"] = x**4 + np.sum(y) + np.sum(z**2)
        funcs["con3"] = x + np.sum(z)

        return funcs, False

    def sens(self, xdict, funcs):
        x = xdict["x"]
        y = xdict["y"]
        z = xdict["z"]

        funcsSens = {
            "obj": {
                "x": 2 * x,
                "y": 4 * y,
                "z": 3 * np.ones(2 * self.N),
            },
            "con1": {
                "x": 2.05 * x * (x * x) ** 0.025,
            },
            "con2": {
                "x": 4 * x**3,
                "y": np.ones(self.N),
                "z": 2 * z,
            },
            "con3": {
                "x": 1.0,
                "z": np.ones(2 * self.N),
            },
        }

        return funcsSens, False

    def setup_optProb(self, sparse=True):
        # set N
        if sparse:
            self.N = 10000
        else:
            self.N = 500

        # Optimization Object
        self.optProb = Optimization("large and sparse", self.objfunc, sens=self.sens)

        # Design Variables
        self.optProb.addVar("x", lower=-100, upper=150, value=0)
        self.optProb.addVarGroup("y", self.N, lower=-10 - np.arange(self.N), upper=np.arange(self.N), value=0)
        self.optProb.addVarGroup(
            "z", 2 * self.N, upper=np.arange(2 * self.N), lower=-100 - np.arange(2 * self.N), value=0
        )

        # Constraints
        self.optProb.addCon("con1", upper=100, wrt=["x"])
        self.optProb.addCon("con2", upper=100)
        self.optProb.addCon("con3", lower=4, wrt=["x", "z"])
        xJac = np.ones((self.N, 1))
        if sparse:
            rows_cols = np.array([i for i in range(0, self.N)]).astype(int)
            yJac = {"coo": [rows_cols, rows_cols, np.ones(self.N)], "shape": [self.N, self.N]}
        else:
            yJac = np.eye(self.N)
        self.optProb.addConGroup(
            "lincon",
            self.N,
            lower=2 - 3 * np.arange(self.N),
            linear=True,
            wrt=["x", "y"],
            jac={"x": xJac, "y": yJac},
        )
        self.optProb.addObj("obj")

    @parameterized.expand(
        [
            ("SNOPT", True),
            ("IPOPT", True),
            ("SNOPT", False),
        ]
    )
    def test_opt(self, optName, sparse):
        self.optName = optName
        self.setup_optProb(sparse=sparse)
        sol = self.optimize()
        self.assert_solution_allclose(sol, 1e-5, partial_x=True)

    def test_dense_workspace_too_small(self):
        self.optName = "SNOPT"
        self.setup_optProb(sparse=False)
        optOptions = {"Total real workspace": 401300}  # 500 + 200 * (503 + 1501)
        sol = self.optimize(optOptions=optOptions)

        # Check that the workspace is too small without overwriting the lengths
        self.assert_inform_equal(sol, 84)


if __name__ == "__main__":
    unittest.main()
