"""Test NSGA2."""

# Standard Python modules
import unittest

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
    histFileName = None

    def objfunc(self, xdict):
        x = xdict["x"]
        y = xdict["y"]

        funcs = {}
        # Adding an offset so that fStar != 0.0
        funcs["obj1"] = (x - 0.0) ** 2 + (y - 0.0) ** 2 + 10
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

    @parameterized.expand(
        [
            (1,),
            (2,),
        ]
    )
    def test_opt(self, n_obj):
        self.setup_optProb(n_obj)

        # 300 generations will find x=(0,0), 200 or less will find x=(1,1)
        optOptions = {"maxGen": 200}

        sol = self.optimize(optOptions=optOptions)
        tol = 1e-2
        if n_obj == 1:
            assert_allclose(sol.xStar["x"], 0.0, atol=tol, rtol=tol)
            assert_allclose(sol.xStar["y"], 0.0, atol=tol, rtol=tol)
            assert_allclose(sol.fStar, 10.0, atol=tol, rtol=tol)
        elif n_obj == 2:
            assert_allclose(sol.xStar["x"], 1.0, atol=tol, rtol=tol)
            assert_allclose(sol.xStar["y"], 1.0, atol=tol, rtol=tol)
            assert_allclose(sol.fStar["obj1"], 12.0, atol=tol, rtol=tol)
            assert_allclose(sol.fStar["obj2"], 0.0, atol=tol, rtol=tol)

    def test_options(self):
        with self.assertRaises(ValueError):
            # PopSize must be a multiple of 4
            self.setup_optProb(1)
            optOptions = {"PopSize": 5}
            self.optimize(optOptions=optOptions)


if __name__ == "__main__":
    unittest.main()
