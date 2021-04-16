"""
This example is taken from the OpenOpt Examples website.

https://github.com/PythonCharmers/OOSuite/blob/master/FuncDesigner/FuncDesigner/examples/nlpSparse.py

Only testing with SNOPT, which supports sparse formats.
"""

# Standard Python modules
import unittest

# External modules
import numpy as np
import scipy

# First party modules
from pyoptsparse import Optimization
from testing_utils import OptTest


class TestLarge(OptTest):
    xStar = None
    fStar = 10.0

    def objfunc(self, xdict):
        x = xdict["x"]
        y = xdict["y"]
        z = xdict["z"]
        funcs = {}
        funcs["obj"] = x ** 2 + 2 * np.sum(y ** 2) + 3 * np.sum(z)
        funcs["con1"] = x + 1e-3 * abs(x) ** 2.05
        funcs["con2"] = x ** 4 + np.sum(y) + np.sum(z ** 2)
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
                "x": 4 * x ** 3,
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
            self.N = 50000
        else:
            self.N = 500
        self.xStar = {"x": 2}
        # Optimization Object
        optProb = Optimization("large and sparse", self.objfunc, sens=self.sens)

        # Design Variables
        optProb.addVar("x", lower=-100, upper=150, value=0)
        optProb.addVarGroup("y", self.N, lower=-10 - np.arange(self.N), upper=np.arange(self.N), value=0)
        optProb.addVarGroup("z", 2 * self.N, upper=np.arange(2 * self.N), lower=-100 - np.arange(2 * self.N), value=0)

        # Constraints
        optProb.addCon("con1", upper=100, wrt=["x"])
        optProb.addCon("con2", upper=100)
        optProb.addCon("con3", lower=4, wrt=["x", "z"])
        xJac = np.ones((self.N, 1))
        if sparse:
            yJac = scipy.sparse.spdiags(np.ones(self.N), 0, self.N, self.N)
        else:
            yJac = np.eye(self.N)
        optProb.addConGroup(
            "lincon",
            self.N,
            lower=2 - 3 * np.arange(self.N),
            linear=True,
            wrt=["x", "y"],
            jac={"x": xJac, "y": yJac},
        )
        optProb.addObj("obj")
        return optProb

    def test_sparse(self):
        self.optName = "SNOPT"
        optProb = self.setup_optProb(sparse=True)
        sol = self.optimize(optProb)
        self.assert_solution(sol, 1e-5, partial=True)

    def test_sparse_IPOPT(self):
        self.optName = "IPOPT"
        optProb = self.setup_optProb(sparse=True)
        sol = self.optimize(optProb)
        self.assert_solution(sol, 1e-5, partial=True)

    def test_dense_default(self):
        self.optName = "SNOPT"
        optProb = self.setup_optProb(sparse=False)
        sol = self.optimize(optProb)
        self.assert_solution(sol, 1e-5, partial=True)

    def test_dense_user(self):
        self.optName = "SNOPT"
        optProb = self.setup_optProb(sparse=False)
        optOptions = {"Total real workspace": 401300}  # 500 + 200 * (503 + 1501)
        sol = self.optimize(optProb, optOptions=optOptions)

        # Check that the workspace is too small without overwriting the lengths
        self.assert_inform(sol, 84)


if __name__ == "__main__":
    unittest.main()
