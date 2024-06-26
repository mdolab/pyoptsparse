"""
Tests that pyoptsparse allows an objective or gradient function to return a "2"
as its fail status.  This status is returned to the Optimizer indicating that
termination has been requested by the user.

The proper response of the pyIPOPT and pySNOPT optimizers are tested.
"""
# Standard Python modules
import unittest

# External modules
import numpy as np
from parameterized import parameterized

# First party modules
from pyoptsparse import OPT, Optimization


class TerminateComp:
    def __init__(self, max_obj=1000, max_sens=1000):
        self.max_obj = max_obj
        self.max_sens = max_sens
        self.obj_count = 0
        self.sens_count = 0

    def objfunc(self, xdict):
        """Evaluates the equation f(x,y) = (x-3)^2 + xy + (y+4)^2 - 3"""
        x = xdict["x"]
        y = xdict["y"]
        funcs = {}

        funcs["obj"] = (x - 3.0) ** 2 + x * y + (y + 4.0) ** 2 - 3.0

        if self.obj_count > self.max_obj:
            fail = 2
        else:
            self.obj_count += 1
            fail = False

        return funcs, fail

    def sens(self, xdict, funcs):
        """f(x,y) = (x-3)^2 + xy + (y+4)^2 - 3"""
        x = xdict["x"]
        y = xdict["y"]

        funcsSens = {
            "obj": {
                "x": 2.0 * x - 6.0 + y,
                "y": 2.0 * y + 8.0 + x,
            }
        }

        if self.sens_count > self.max_sens:
            fail = 2
        else:
            self.sens_count += 1
            fail = False

        return funcsSens, fail


def setup_optProb(termcomp):
    optProb = Optimization("Paraboloid", termcomp.objfunc)

    optProb.addVarGroup("x", 1, varType="c", lower=-50.0, upper=50.0, value=0.0)
    optProb.addVarGroup("y", 1, varType="c", lower=-50.0, upper=50.0, value=0.0)

    optProb.addObj("obj")

    con_jac = {"x": np.array(-1.0), "y": np.array(1.0)}

    optProb.addConGroup("con", 1, lower=-15.0, upper=-15.0, wrt=["x", "y"], linear=True, jac=con_jac)

    return optProb


class TestUserTerminationStatus(unittest.TestCase):
    optOptions = {
        "IPOPT": {
            "output_file": "{}.out",
        },
        "SNOPT": {
            "Print file": "{}.out",
            "Summary file": "{}_summary.out",
        },
    }

    optExitCode = {
        "IPOPT": 5,
        "SNOPT": 71,
    }

    @parameterized.expand(["IPOPT", "SNOPT"])
    def test_obj(self, optName):
        termcomp = TerminateComp(max_obj=2)
        optProb = setup_optProb(termcomp)

        optOptions = self.optOptions[optName]
        for key, val in optOptions.items():
            optOptions[key] = val.format(f"{optName}_user_terminate_obj")

        try:
            opt = OPT(optName, options=optOptions)
        except ImportError:
            raise unittest.SkipTest(f"Optimizer not available: {optName}")

        sol = opt(optProb, sens=termcomp.sens)

        self.assertEqual(termcomp.obj_count, 3)

        # Exit code for user requested termination.
        self.assertEqual(sol.optInform["value"], self.optExitCode[optName])

    @parameterized.expand(["IPOPT", "SNOPT"])
    def test_sens(self, optName):
        termcomp = TerminateComp(max_sens=3)
        optProb = setup_optProb(termcomp)

        optOptions = self.optOptions[optName]
        for key, val in optOptions.items():
            optOptions[key] = val.format(f"{optName}_user_terminate_sens")

        try:
            opt = OPT(optName, options=optOptions)
        except ImportError:
            raise unittest.SkipTest(f"Optimizer not available: {optName}")

        sol = opt(optProb, sens=termcomp.sens)

        self.assertEqual(termcomp.sens_count, 4)

        # Exit code for user requested termination.
        self.assertEqual(sol.optInform["value"], self.optExitCode[optName])


if __name__ == "__main__":
    unittest.main()
