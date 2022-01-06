""" Test NSGA2."""

# Standard Python modules
import unittest

# External modules
from numpy.testing import assert_allclose

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

    def setup_optProb(self):
        # Instantiate Optimization Problem
        self.optProb = Optimization("quadratic", self.objfunc)
        self.optProb.addVar("x", value=0, lower=-600, upper=600)
        self.optProb.addVar("y", value=0, lower=-600, upper=600)

        self.optProb.addObj("obj1")
        self.optProb.addObj("obj2")

    def test_opt(self):
        self.setup_optProb()

        # 300 generations will find x=(0,0), 200 or less will find x=(1,1)
        optOptions = {"maxGen": 200}
        sol = self.optimize(optOptions=optOptions)
        tol = 1e-2
        assert_allclose(sol.variables["x"][0].value, 1.0, atol=tol, rtol=tol)
        assert_allclose(sol.variables["y"][0].value, 1.0, atol=tol, rtol=tol)


if __name__ == "__main__":
    unittest.main()
