"""
Regression test for a DV ``offset`` combined with a linear constraint.

An active linear constraint must be enforced about the correct intercept regardless of a DV
offset. This exercises both SNOPT, which evaluates linear-constraint rows internally as
``jac @ x_opt``, and the other wrappers, which evaluate linear constraints in user space via
``evaluateLinearConstraints``.
"""

# Standard Python modules
import unittest

# External modules
import numpy as np
from numpy.testing import assert_allclose
from parameterized import parameterized

# First party modules
from pyoptsparse import OPT, Optimization


def objfunc(xdict):
    """min (x-2)^2 + (y-2)^2 with an inactive nonlinear constraint x^2 + y^2 <= 100."""
    x = xdict["xvars"]
    funcs = {}
    funcs["obj"] = (x[0] - 2.0) ** 2 + (x[1] - 2.0) ** 2
    funcs["nlcon"] = np.array([x[0] ** 2 + x[1] ** 2])
    return funcs, False


def sens(xdict, funcs):
    x = xdict["xvars"]
    funcsSens = {
        "obj": {"xvars": np.array([2.0 * (x[0] - 2.0), 2.0 * (x[1] - 2.0)])},
        "nlcon": {"xvars": np.array([[2.0 * x[0], 2.0 * x[1]]])},
    }
    return funcsSens, False


class TestLinearConstraintOffset(unittest.TestCase):
    # The linear constraint x + y <= 1 is active at the true optimum (0.5, 0.5), so a mishandled
    # DV offset shifts the enforced bound and moves the solution. The inactive nonlinear constraint
    # is required for SNOPT to route the linear row through its internal evaluation (otherwise the
    # row is evaluated in user space and the offset handling is never exercised).

    @staticmethod
    def options(optName):
        if optName == "SNOPT":
            return {
                "Major feasibility tolerance": 1e-8,
                "Major optimality tolerance": 1e-8,
                "Print file": "lincon_offset_SNOPT.out",
                "Summary file": "lincon_offset_SNOPT_summary.out",
            }
        if optName == "IPOPT":
            return {"print_level": 0, "output_file": "lincon_offset_IPOPT.out"}
        return {}

    def optimize(self, optName, offset):
        optProb = Optimization("lincon_offset", objfunc)
        optProb.addVarGroup("xvars", 2, lower=-50.0, upper=50.0, value=0.0, offset=offset)
        optProb.addObj("obj")
        optProb.addConGroup("nlcon", 1, upper=100.0, wrt=["xvars"])
        optProb.addConGroup(
            "lincon", 1, upper=1.0, wrt=["xvars"], linear=True, jac={"xvars": np.array([[1.0, 1.0]])}
        )
        try:
            opt = OPT(optName, options=self.options(optName))
        except ImportError as e:
            raise unittest.SkipTest(f"Optimizer not available: {optName}") from e
        return opt(optProb, sens=sens)

    @parameterized.expand(["SNOPT", "SLSQP", "IPOPT"])
    def test_offset_invariance(self, optName):
        # The optimum must be the same with or without a DV offset.
        assert_allclose(self.optimize(optName, 0.0).xStar["xvars"], [0.5, 0.5], atol=1e-5, rtol=1e-5)
        assert_allclose(self.optimize(optName, 3.0).xStar["xvars"], [0.5, 0.5], atol=1e-5, rtol=1e-5)


if __name__ == "__main__":
    unittest.main()
