"""Test class for SLSQP specific tests"""

# Standard Python modules
import unittest

# First party modules
from pyoptsparse import OPT, Optimization


class TestSLSQP(unittest.TestCase):
    def test_slsqp_strong_bound_enforcement(self):
        """
        Test that SLSQP will never evaluate the function or gradient outside
        the design variable bounds. Without strong bound enforcement, the
        optimizer will step outside the bounds and a ValueError will be raised.
        With strong bound enforement, this code will run without raising any
        errors
        """

        def objfunc(xdict):
            x = xdict["xvars"]
            funcs = {}
            if x[0] < 0:
                raise ValueError("Function cannot be evaluated below 0.")
            funcs["obj"] = (x[0] + 1.0) ** 2
            fail = False
            return funcs, fail

        def sens(xdict, funcs):
            x = xdict["xvars"]
            if x[0] < 0:
                raise ValueError("Function cannot be evaluated below 0.")
            funcsSens = {
                "obj": {"xvars": [2 * (x[0] + 1.0)]},
            }
            fail = False
            return funcsSens, fail

        optProb = Optimization("Problem with Error Region", objfunc)
        optProb.addVarGroup("xvars", 1, lower=[0], value=[2])
        optProb.addObj("obj")
        opt = OPT("SLSQP")
        sol = opt(optProb, sens=sens)
        self.assertEqual(sol.optInform["value"], 0)
        self.assertGreaterEqual(sol.xStar["xvars"][0], 0)
        self.assertAlmostEqual(sol.xStar["xvars"][0], 0, places=9)
