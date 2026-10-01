"""
Regression tests for a DV ``offset`` combined with linear constraints.

An active linear constraint must be enforced about the correct intercept regardless of a DV offset.
SNOPT is the only backend that evaluates linear-constraint rows internally (as ``A @ x_opt`` in
optimizer space, where a DV offset shifts the row value); the other wrappers evaluate linear
constraints in user space via ``evaluateLinearConstraints`` (where the offset is already applied).
These tests cover both code paths, and both of SNOPT's sub-cases:

* ``test_offset_with_nonlinear_constraint``: a nonlinear constraint is present, so SNOPT keeps the
  linear row on its internal ``A`` path. This is the primary trigger for the bug.
* ``test_offset_all_linear``: an all-linear problem. SNOPT cannot have zero nonlinear constraints,
  so it reclassifies the *first* constraint as a dummy nonlinear one (evaluated in user space via the
  callback) while still evaluating the *remaining* linear rows from ``A``. The bug therefore survives
  on every linear constraint except the first, so this needs at least two linear constraints to
  exercise it.

In every case the optimum must be invariant to the offset.
"""

# Standard Python modules
import unittest

# External modules
import numpy as np
from numpy.testing import assert_allclose
from parameterized import parameterized

# First party modules
from pyoptsparse import OPT, Optimization


def _objective(x):
    """Paraboloid centered at (2, 2); returns the value and its gradient."""
    value = (x[0] - 2.0) ** 2 + (x[1] - 2.0) ** 2
    grad = np.array([2.0 * (x[0] - 2.0), 2.0 * (x[1] - 2.0)])
    return value, grad


def nonlinear_objfunc(xdict):
    x = xdict["x"]
    value, _ = _objective(x)
    return {"obj": value, "nlcon": np.array([x[0] ** 2 + x[1] ** 2])}, False


def nonlinear_sens(xdict, funcs):
    x = xdict["x"]
    _, grad = _objective(x)
    return {"obj": {"x": grad}, "nlcon": {"x": np.array([[2.0 * x[0], 2.0 * x[1]]])}}, False


def linear_objfunc(xdict):
    value, _ = _objective(xdict["x"])
    return {"obj": value}, False


def linear_sens(xdict, funcs):
    _, grad = _objective(xdict["x"])
    return {"obj": {"x": grad}}, False


def options(optName, tag):
    if optName == "SNOPT":
        return {
            "Major feasibility tolerance": 1e-8,
            "Major optimality tolerance": 1e-8,
            "Print file": f"lincon_offset_{tag}_SNOPT.out",
            "Summary file": f"lincon_offset_{tag}_SNOPT_summary.out",
        }
    if optName == "IPOPT":
        return {"print_level": 0, "output_file": f"lincon_offset_{tag}_IPOPT.out"}
    return {}


class TestLinearConstraintOffset(unittest.TestCase):
    """A DV offset must not move the optimum of a problem with an active linear constraint.

    Both problems minimize the paraboloid centered at (2, 2) subject to the active linear constraint
    ``x + y <= 1``, whose true optimum is (0.5, 0.5). If the offset shifts the enforced linear bound,
    the constraint effectively deactivates and the solver returns the unconstrained minimum (2, 2).
    """

    @staticmethod
    def optimize(optName, optProb, sens, tag):
        try:
            opt = OPT(optName, options=options(optName, tag))
        except ImportError as e:
            raise unittest.SkipTest(f"Optimizer not available: {optName}") from e
        return opt(optProb, sens=sens)

    @staticmethod
    def nonlinear_plus_linear(offset):
        # The inactive nonlinear constraint keeps SNOPT treating the problem as nonlinear, so the
        # linear row stays on SNOPT's internal A @ x_opt path (where the bug lives). Without it, SNOPT
        # would route the lone linear row through the user-space callback and the bug would be masked.
        optProb = Optimization("nonlinear_plus_linear", nonlinear_objfunc)
        optProb.addVarGroup("x", 2, lower=-50.0, upper=50.0, value=0.0, offset=offset)
        optProb.addObj("obj")
        optProb.addConGroup("nlcon", 1, upper=100.0, wrt=["x"])
        optProb.addConGroup("lincon", 1, upper=1.0, wrt=["x"], linear=True, jac={"x": np.array([[1.0, 1.0]])})
        return optProb

    @staticmethod
    def two_linear(offset):
        # All-linear: SNOPT makes the first constraint (c0) a dummy nonlinear one evaluated in user
        # space, but still evaluates c1 from A. c0 is loose/inactive; the binding c1 (x + y <= 1) must
        # be enforced correctly under the offset, i.e. the offset handling must reach linear rows
        # beyond the dummy.
        optProb = Optimization("two_linear", linear_objfunc)
        optProb.addVarGroup("x", 2, lower=-50.0, upper=50.0, value=0.0, offset=offset)
        optProb.addObj("obj")
        optProb.addConGroup("c0", 1, upper=100.0, wrt=["x"], linear=True, jac={"x": np.array([[1.0, -1.0]])})
        optProb.addConGroup("c1", 1, upper=1.0, wrt=["x"], linear=True, jac={"x": np.array([[1.0, 1.0]])})
        return optProb

    @parameterized.expand(["SNOPT", "SLSQP", "IPOPT"])
    def test_offset_with_nonlinear_constraint(self, optName):
        for offset in (0.0, 3.0):
            sol = self.optimize(optName, self.nonlinear_plus_linear(offset), nonlinear_sens, "nl")
            assert_allclose(sol.xStar["x"], [0.5, 0.5], atol=1e-5, rtol=1e-5)

    @parameterized.expand(["SNOPT", "SLSQP", "IPOPT"])
    def test_offset_all_linear(self, optName):
        for offset in (0.0, 3.0):
            sol = self.optimize(optName, self.two_linear(offset), linear_sens, "lin")
            assert_allclose(sol.xStar["x"], [0.5, 0.5], atol=1e-5, rtol=1e-5)


if __name__ == "__main__":
    unittest.main()
