"""
Tests that pyoptsparse allows a objective or gradient function to return a "2" as its fail status.
This status is passed to pySNOPT, which returns a -2 fail status, indicating that termination has
been requested by the user.
"""
# Standard Python modules
import unittest

# External modules
import numpy as np

# First party modules
from pyoptsparse import SNOPT, Optimization, IPOPT
from pyoptsparse.pyOpt_error import Error


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
        conval = -x + y
        funcs["con"] = conval

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


con_jac = {}
con_jac["x"] = np.array(-1.0)
con_jac["y"] = np.array(1.0)


class TestUserTerminationStatus(unittest.TestCase):

    def test_obj_SNOPT(self):
        termcomp = TerminateComp(max_obj=2)
        optProb = Optimization("Paraboloid", termcomp.objfunc)

        optProb.addVarGroup("x", 1, varType="c", lower=-50.0, upper=50.0, value=0.0)
        optProb.addVarGroup("y", 1, varType="c", lower=-50.0, upper=50.0, value=0.0)

        optProb.addObj("obj")

        optProb.addConGroup("con", 1, lower=-15.0, upper=-15.0, wrt=["x", "y"], linear=True, jac=con_jac)

        test_name = "SNOPT_user_termination_obj"
        optOptions = {
            "Print file": f"{test_name}.out",
            "Summary file": f"{test_name}_summary.out",
        }
        try:
            opt = SNOPT(options=optOptions)
        except Error as e:
            if "There was an error importing" in e.message:
                raise unittest.SkipTest("Optimizer not available: SNOPT")
            raise e

        sol = opt(optProb, sens=termcomp.sens)

        self.assertEqual(termcomp.obj_count, 3)

        # Exit code for user requested termination.
        self.assertEqual(sol.optInform["value"], 71)

    def test_sens_SNOPT(self):
        termcomp = TerminateComp(max_sens=3)
        optProb = Optimization("Paraboloid", termcomp.objfunc)

        optProb.addVarGroup("x", 1, varType="c", lower=-50.0, upper=50.0, value=0.0)
        optProb.addVarGroup("y", 1, varType="c", lower=-50.0, upper=50.0, value=0.0)

        optProb.addObj("obj")

        optProb.addConGroup("con", 1, lower=-15.0, upper=-15.0, wrt=["x", "y"], linear=True, jac=con_jac)

        test_name = "SNOPT_user_termination_sens"
        optOptions = {
            "Print file": f"{test_name}.out",
            "Summary file": f"{test_name}_summary.out",
        }
        try:
            opt = SNOPT(options=optOptions)
        except Error as e:
            if "There was an error importing" in e.message:
                raise unittest.SkipTest("Optimizer not available: SNOPT")
            raise e

        sol = opt(optProb, sens=termcomp.sens)

        self.assertEqual(termcomp.sens_count, 4)

        # Exit code for user requested termination.
        self.assertEqual(sol.optInform["value"], 71)

    def test_obj_IPOPT(self):
        termcomp = TerminateComp(max_obj=2)
        optProb = Optimization("Paraboloid", termcomp.objfunc)

        optProb.addVarGroup("x", 1, varType="c", lower=-50.0, upper=50.0, value=0.0)
        optProb.addVarGroup("y", 1, varType="c", lower=-50.0, upper=50.0, value=0.0)

        optProb.addObj("obj")

        optProb.addConGroup("con", 1, lower=-15.0, upper=-15.0, wrt=["x", "y"], linear=True, jac=con_jac)

        test_name = "IPOPT_user_termination_obj"
        optOptions = {
            "output_file": f"{test_name}.out",
        }
        try:
            opt = IPOPT(options=optOptions)
        except Error as e:
            if "There was an error importing" in e.message:
                raise unittest.SkipTest("Optimizer not available: IPOPT")
            raise e

        sol = opt(optProb, sens=termcomp.sens)

        self.assertEqual(termcomp.obj_count, 3)

        # Exit code for user requested termination.
        self.assertEqual(sol.optInform["value"], 5)

    def test_sens_IPOPT(self):
        termcomp = TerminateComp(max_sens=3)
        optProb = Optimization("Paraboloid", termcomp.objfunc)

        optProb.addVarGroup("x", 1, varType="c", lower=-50.0, upper=50.0, value=0.0)
        optProb.addVarGroup("y", 1, varType="c", lower=-50.0, upper=50.0, value=0.0)

        optProb.addObj("obj")

        optProb.addConGroup("con", 1, lower=-15.0, upper=-15.0, wrt=["x", "y"], linear=True, jac=con_jac)

        test_name = "IPOPT_user_termination_sens"
        optOptions = {
            "output_file": f"{test_name}.out",
        }
        try:
            opt = IPOPT(options=optOptions)
        except Error as e:
            if "There was an error importing" in e.message:
                raise unittest.SkipTest("Optimizer not available: IPOPT")
            raise e

        sol = opt(optProb, sens=termcomp.sens)

        self.assertEqual(termcomp.sens_count, 4)

        # Exit code for user requested termination.
        self.assertEqual(sol.optInform["value"], 5)


if __name__ == "__main__":
    unittest.main()
