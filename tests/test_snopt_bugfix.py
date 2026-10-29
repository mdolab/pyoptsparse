# Regression tests for pySNOPT bugs, all on the paraboloid
# f(x, y) = (x-3)^2 + xy + (y+4)^2 - 3 with no, or only linear, constraints.

# Standard Python modules
import unittest

# External modules
import numpy as np
from numpy.testing import assert_allclose

# First party modules
from pyoptsparse import SNOPT, Optimization


def objfunc(xdict):
    """Evaluates the equation f(x,y) = (x-3)^2 + xy + (y+4)^2 - 3"""
    x = xdict["x"]
    y = xdict["y"]
    funcs = {}

    funcs["obj"] = (x - 3.0) ** 2 + x * y + (y + 4.0) ** 2 - 3.0

    fail = False
    return funcs, fail


def sens(xdict, funcs):
    """f(x,y) = (x-3)^2 + xy + (y+4)^2 - 3"""
    x = xdict["x"]
    y = xdict["y"]

    funcsSens = {
        "obj": {
            "x": 2.0 * x - 6.0 + y,
            "y": 2.0 * y + 8.0 + x,
        }
    }

    fail = False
    return funcsSens, fail


con_jac = {}
con_jac["x"] = np.array(-1.0)
con_jac["y"] = np.array(1.0)


class TestSNOPTBug(unittest.TestCase):
    def test_opt(self):
        # pySNOPT wouldn't optimize a model that has only equality constraints because it thought the problem was
        # trivial. The minimum should be at (7.166667, -7.833334), but with the bug, x and y stay at zero.
        # Optimization Object
        optProb = Optimization("Paraboloid", objfunc)

        # Design Variables
        optProb.addVarGroup("x", 1, varType="c", lower=-50.0, upper=50.0, value=0.0)
        optProb.addVarGroup("y", 1, varType="c", lower=-50.0, upper=50.0, value=0.0)

        # Objective
        optProb.addObj("obj")

        # Equality Constraint
        optProb.addConGroup("con", 1, lower=-15.0, upper=-15.0, wrt=["x", "y"], linear=True, jac=con_jac)

        # Check optimization problem:
        print(optProb)
        test_name = "bugfix_SNOPT_test_opt"
        optOptions = {
            "Major feasibility tolerance": 1e-1,
            "Print file": f"{test_name}.out",
            "Summary file": f"{test_name}_summary.out",
        }

        # Optimizer
        try:
            opt = SNOPT(options=optOptions)
        except ImportError as e:
            raise unittest.SkipTest("Optimizer not available: SNOPT") from e

        sol = opt(optProb, sens=sens)

        # Check Solution 7.166667, -7.833334
        tol = 1e-6
        assert_allclose(sol.xStar["x"], 7.166667, atol=tol, rtol=tol)
        assert_allclose(sol.xStar["y"], -7.833333, atol=tol, rtol=tol)

    def test_opt_bug1(self):
        # Due to a new feature, there is a TypeError when you optimize a model without a constraint.
        optProb = Optimization("Paraboloid", objfunc)

        # Design Variables
        optProb.addVarGroup("x", 1, varType="c", lower=-50.0, upper=50.0, value=0.0)
        optProb.addVarGroup("y", 1, varType="c", lower=-50.0, upper=50.0, value=0.0)

        # Objective
        optProb.addObj("obj")

        test_name = "bugfix_SNOPT_bug1"
        optOptions = {
            "Major feasibility tolerance": 1e-1,
            "Print file": f"{test_name}.out",
            "Summary file": f"{test_name}_summary.out",
        }

        # Optimizer
        try:
            opt = SNOPT(options=optOptions)
        except ImportError as e:
            raise unittest.SkipTest("Optimizer not available: SNOPT") from e

        sol = opt(optProb, sens=sens)

        # Check unconstrained minimum (20/3, -22/3)
        tol = 1e-6
        assert_allclose(sol.xStar["x"], 20.0 / 3.0, atol=tol, rtol=tol)
        assert_allclose(sol.xStar["y"], -22.0 / 3.0, atol=tol, rtol=tol)

    def test_opt_bug_print_2con(self):
        # Optimization Object
        optProb = Optimization("Paraboloid", objfunc)

        # Design Variables
        optProb.addVarGroup("x", 1, varType="c", lower=-50.0, upper=50.0, value=0.0)
        optProb.addVarGroup("y", 1, varType="c", lower=-50.0, upper=50.0, value=0.0)

        # Objective
        optProb.addObj("obj")

        con_jac2 = {}
        con_jac2["x"] = -np.ones((2, 1))
        con_jac2["y"] = np.ones((2, 1))

        con_jac3 = {}
        con_jac3["x"] = -np.ones((3, 1))
        con_jac3["y"] = np.ones((3, 1))

        # Equality Constraint
        optProb.addConGroup("con", 2, lower=-15.0, upper=-15.0, wrt=["x", "y"], linear=True, jac=con_jac2)
        optProb.addConGroup("con2", 3, lower=-15.0, upper=-15.0, wrt=["x", "y"], linear=True, jac=con_jac3)

        # Check optimization problem:
        print(optProb)

        test_name = "bugfix_SNOPT_bug_print_2con"
        optOptions = {
            "Major feasibility tolerance": 1e-1,
            "Print file": f"{test_name}.out",
            "Summary file": f"{test_name}_summary.out",
        }

        # Optimizer
        try:
            opt = SNOPT(options=optOptions)
        except ImportError as e:
            raise unittest.SkipTest("Optimizer not available: SNOPT") from e

        sol = opt(optProb, sens=sens)

        print(sol)

        # Every row is -x + y = -15, the same constraint as test_opt, so the solution is the same
        tol = 1e-6
        assert_allclose(sol.xStar["x"], 7.166667, atol=tol, rtol=tol)
        assert_allclose(sol.xStar["y"], -7.833333, atol=tol, rtol=tol)

    def test_opt_all_linear_inequality_before_equality(self):
        # For an all-linear problem SNOPT makes its first row a dummy nonlinear constraint evaluated via the callback.
        # That row is the first in the ["ne", "ni", "le", "li"] ordering (the equality), not the first added (the
        # inequality). Pairing the inequality value with the equality bounds enforces x = 0 and drops x = y, giving
        # (0, -4) instead of the true optimum (0, 0).
        optProb = Optimization("Paraboloid", objfunc)

        # Design Variables
        optProb.addVarGroup("x", 1, varType="c", lower=-50.0, upper=50.0, value=0.0)
        optProb.addVarGroup("y", 1, varType="c", lower=-50.0, upper=50.0, value=0.0)

        # Objective
        optProb.addObj("obj")

        # Inequality x >= 0 added before equality x - y = 0
        optProb.addConGroup("ineq", 1, lower=0.0, wrt=["x"], linear=True, jac={"x": np.array([[1.0]])})
        optProb.addConGroup("eq", 1, lower=0.0, upper=0.0, wrt=["x", "y"], linear=True, jac={"x": 1.0, "y": -1.0})

        test_name = "bugfix_SNOPT_all_linear_inequality_before_equality"
        optOptions = {
            "Major feasibility tolerance": 1e-8,
            "Major optimality tolerance": 1e-8,
            "Print file": f"{test_name}.out",
            "Summary file": f"{test_name}_summary.out",
        }

        # Optimizer
        try:
            opt = SNOPT(options=optOptions)
        except ImportError as e:
            raise unittest.SkipTest("Optimizer not available: SNOPT") from e

        sol = opt(optProb, sens=sens)

        tol = 1e-6
        assert_allclose(sol.xStar["x"], 0.0, atol=tol, rtol=tol)
        assert_allclose(sol.xStar["y"], 0.0, atol=tol, rtol=tol)
        assert_allclose(sol.fStar, 22.0, atol=tol, rtol=tol)


if __name__ == "__main__":
    unittest.main()
