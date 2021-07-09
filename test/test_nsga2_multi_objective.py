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


if __name__ == "__main__":
    unittest.main()
