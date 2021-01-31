"""
This example is taken from the OpenOpt Examples website.

https://github.com/PythonCharmers/OOSuite/blob/master/FuncDesigner/FuncDesigner/examples/nlpSparse.py

Only testing with SNOPT, which supports sparse formats.
"""

import unittest
import numpy as np
from numpy.testing import assert_allclose
from numpy import arange
import scipy

from pyoptsparse import Optimization, SNOPT
from pyoptsparse.pyOpt_error import Error


class TestLarge(unittest.TestCase):
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

    def optimize(self, sparse=True, tol=None, optOptions={}, storeHistory=False):
        # set N
        if sparse:
            self.N = 50000
        else:
            self.N = 500
        # Optimization Object
        optProb = Optimization("large and sparse", self.objfunc)

        # Design Variables
        optProb.addVar("x", lower=-100, upper=150, value=0)
        optProb.addVarGroup("y", self.N, lower=-10 - arange(self.N), upper=arange(self.N), value=0)
        optProb.addVarGroup("z", 2 * self.N, upper=arange(2 * self.N), lower=-100 - arange(2 * self.N), value=0)

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
            lower=2 - 3 * arange(self.N),
            linear=True,
            wrt=["x", "y"],
            jac={"x": xJac, "y": yJac},
        )
        optProb.addObj("obj")

        # Optimizer
        try:
            opt = SNOPT(options=optOptions)
        except Error:
            raise unittest.SkipTest("Optimizer not available: SNOPT")

        sol = opt(optProb, sens=self.sens)

        # Check Solution
        if tol is not None:
            assert_allclose(sol.objectives["obj"].value, 10.0, atol=tol, rtol=tol)
            assert_allclose(sol.variables["x"][0].value, 2.0, atol=tol, rtol=tol)
        return sol

    def test_sparse(self):
        test_name = "large_sparse_SNOPT"
        optOptions = {
            "Print file": "{}.out".format(test_name),
            "Summary file": "{}_summary.out".format(test_name),
        }
        self.optimize(tol=1e-5, optOptions=optOptions)

    def test_default(self):
        test_name = "SNOPT_workspace_default"
        optOptions = {
            "Print file": "{}.out".format(test_name),
            "Summary file": "{}_summary.out".format(test_name),
        }
        sol = self.optimize(sparse=False, optOptions=optOptions)

        # Check that overwriting the workspace lengths works
        tol = 1e-5
        assert_allclose(sol.objectives["obj"].value, 10.0, atol=tol, rtol=tol)
        assert_allclose(sol.variables["x"][0].value, 2.0, atol=tol, rtol=tol)

    def test_user(self):
        test_name = "SNOPT_workspace_user"
        optOptions = {
            "Print file": "{}.out".format(test_name),
            "Summary file": "{}_summary.out".format(test_name),
            "Total real workspace": 401300,  # 500 + 200 * (503 + 1501)
        }
        sol = self.optimize(sparse=False, optOptions=optOptions)

        # Check that the workspace is too small without overwriting the lengths
        self.assertEqual(sol.optInform["value"], 84)


if __name__ == "__main__":
    unittest.main()
