""" Test NSGA2."""

# Standard Python modules
import sys
import unittest
import warnings

# External modules
from numpy.testing import assert_allclose
from parameterized import parameterized

# First party modules
from pyoptsparse import Optimization

# Local modules
from testing_utils import OptTest


class TestNSGA2(OptTest):
    name = "quadratic"
    optName = "NSGA2"

    def objfunc(self, xdict):
        x = xdict["x"]
        y = xdict["y"]

        funcs = {}
        funcs["obj1"] = (x - 0.0) ** 2 + (y - 0.0) ** 2
        funcs["obj2"] = (x - 1.0) ** 2 + (y - 1.0) ** 2

        fail = False

        return funcs, fail

    def setup_optProb(self, n_obj):
        # Instantiate Optimization Problem
        self.optProb = Optimization("quadratic", self.objfunc)
        self.optProb.addVar("x", value=0, lower=-600, upper=600)
        self.optProb.addVar("y", value=0, lower=-600, upper=600)

        self.optProb.addObj("obj1")
        if n_obj == 2:
            self.optProb.addObj("obj2")

    @parameterized.expand([(1,), (2,)])
    def test_opt(self, n_obj):
        self.setup_optProb(n_obj)

        # 300 generations will find x=(0,0), 200 or less will find x=(1,1)
        optOptions = {"maxGen": 200}
        if sys.platform == "win32":
            warnings.warn(
                "test_nsga2_multi_objective.py fails on windows with two objectives! Skipping for now.", stacklevel=2
            )
            return
        sol = self.optimize(optOptions=optOptions)
        tol = 1e-2
        if n_obj == 1:
            assert_allclose(sol.xStar["x"], 0.0, atol=tol, rtol=tol)
            assert_allclose(sol.xStar["y"], 0.0, atol=tol, rtol=tol)
            assert_allclose(sol.fStar, 0.0, atol=tol, rtol=tol)
        elif n_obj == 2:
            assert_allclose(sol.xStar["x"], 1.0, atol=tol, rtol=tol)
            assert_allclose(sol.xStar["y"], 1.0, atol=tol, rtol=tol)
            assert_allclose(sol.fStar["obj1"], 2.0, atol=tol, rtol=tol)
            assert_allclose(sol.fStar["obj2"], 0.0, atol=tol, rtol=tol)


if __name__ == "__main__":
    unittest.main()
