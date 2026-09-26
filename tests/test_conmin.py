"""Test class for CONMIN specific tests"""

# Standard Python modules
import unittest

# External modules
import numpy as np
from numpy.testing import assert_allclose

# First party modules
from pyoptsparse import OPT, Optimization


class TestCONMIN(unittest.TestCase):
    @staticmethod
    def test_nonlinear_equality():
        """CONMIN has no native equality handling, so equality constraints are split into two one-sided inequalities.
        Here we construct a problem with an equality constraint to exercise CONMIN.
        """

        def objfunc(xdict):
            x = xdict["x"]
            funcs = {"obj": (x[0] - 1) ** 2 + (x[1] - 1) ** 2, "con": [x[0] ** 2 + x[1] ** 2]}
            return funcs, False

        optProb = Optimization("CONMIN nonlinear equality", objfunc)
        optProb.addVarGroup("x", 2, value=[1.0, 1.0], lower=-20, upper=20)
        optProb.addConGroup("con", 1, lower=100.0, upper=100.0)
        optProb.addObj("obj")

        opt = OPT("CONMIN", options={"DELFUN": 1e-10, "DABFUN": 1e-10})
        sol = opt(optProb, sens="FD")

        # nearest point on the circle to (1, 1): (1, 1) / sqrt(2) * 10
        assert_allclose(sol.xStar["x"], [5 * np.sqrt(2), 5 * np.sqrt(2)], atol=1e-4)


if __name__ == "__main__":
    unittest.main()
