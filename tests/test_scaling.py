"""
Unit tests for the user <-> optimizer scaling/mapping layer in
pyOpt_optimization.

The existing test_optProb.py exercises the *value* mappings for design
variables, objectives, and constraints. It does NOT touch the gradient and
Jacobian mappings (_mapObjGradtoOpt / _mapConJactoOpt), which combine row
scaling (by the objective/constraint scale) with column scaling (by invXScale,
the chain-rule factor for the DV change of variables). CLAUDE.md flags this as
the single most bug-prone seam in the codebase, so we pin it down directly here.

We only need finalize() (not a full optimization run) since that is what
populates invXScale, conScale, xOffset, and objectiveIdx.
"""

# Standard Python modules
import unittest

# External modules
import numpy as np
from numpy.testing import assert_allclose

# First party modules
from pyoptsparse import Optimization
from pyoptsparse.pyOpt_utils import convertToCSR, convertToDense


class TestScalingMaps(unittest.TestCase):
    def setUp(self):
        # Distinct, non-trivial per-element scales and offsets so that any
        # mixed-up indexing or row/column confusion shows up.
        self.xScale = {"x": np.array([2.0, 0.5, 4.0]), "y": np.array([10.0, 0.1])}
        self.xOffset = {"x": np.array([1.0, -2.0, 0.5]), "y": np.array([0.0, 3.0])}
        self.objScale = 3.0
        self.conScaleVals = {"c1": np.array([5.0, 0.2]), "c2": np.array([7.0])}

        def objfunc(xdict):
            # Never actually called in these tests, but required by the API.
            return {"obj": 0.0, "c1": np.zeros(2), "c2": np.zeros(1)}, False

        optProb = Optimization("scaling-test", objfunc)
        optProb.addVarGroup("x", 3, lower=-10, upper=10, scale=self.xScale["x"], offset=self.xOffset["x"])
        optProb.addVarGroup("y", 2, lower=-10, upper=10, scale=self.xScale["y"], offset=self.xOffset["y"])
        optProb.addObj("obj", scale=self.objScale)
        optProb.addConGroup("c1", 2, lower=-1, upper=1, scale=self.conScaleVals["c1"])
        optProb.addConGroup("c2", 1, lower=-1, upper=1, scale=self.conScaleVals["c2"])
        optProb.finalize()

        self.optProb = optProb
        self.ndvs = optProb.ndvs
        self.nCon = optProb.nCon
        # invXScale = 1/scale, in DV order (x then y)
        self.invXScale = optProb.invXScale
        # conScale in natural (un-reordered) order: c1, c2
        self.conScale = optProb.conScale

    def test_finalize_populated_scales(self):
        assert_allclose(self.invXScale, 1.0 / np.array([2.0, 0.5, 4.0, 10.0, 0.1]))
        assert_allclose(self.conScale, np.array([5.0, 0.2, 7.0]))
        assert_allclose(self.optProb.xOffset, np.array([1.0, -2.0, 0.5, 0.0, 3.0]))

    def test_mapX_roundtrip_and_formula(self):
        rng = np.random.default_rng(0)
        x_user = rng.uniform(-5, 5, self.ndvs)
        x_opt = self.optProb._mapXtoOpt(x_user)
        # x_opt = (x_user - offset) / invXScale
        assert_allclose(x_opt, (x_user - self.optProb.xOffset) / self.invXScale)
        # round trip
        assert_allclose(self.optProb._mapXtoUser(x_opt), x_user)

    def test_mapObjGrad(self):
        # Objective gradient mapping: g_opt = g_user * s_f * invXScale (column/chain-rule scaling).
        rng = np.random.default_rng(1)
        gobj = rng.uniform(-3, 3, (self.optProb.nObj, self.ndvs))
        gobj_orig = gobj.copy()
        gobj_opt = self.optProb._mapObjGradtoOpt(gobj)
        assert_allclose(gobj_opt, gobj * self.objScale * self.invXScale)
        # the method must not mutate its input
        assert_allclose(gobj, gobj_orig)

    def test_mapConJac_formula_and_roundtrip(self):
        # Build an arbitrary dense Jacobian of the right shape and convert to CSR.
        rng = np.random.default_rng(2)
        dense = rng.uniform(-2, 2, (self.nCon, self.ndvs))
        jac = convertToCSR(dense)

        # _mapConJactoOpt works in place: J_opt = diag(conScale) . J . diag(invXScale)
        self.optProb._mapConJactoOpt(jac)
        expected = np.diag(self.conScale) @ dense @ np.diag(self.invXScale)
        assert_allclose(convertToDense(jac), expected)

        # _mapConJactoUser must invert it back to the original.
        self.optProb._mapConJactoUser(jac)
        assert_allclose(convertToDense(jac), dense)

    def test_mapObj_value_roundtrip(self):
        f_user = np.array([2.5])
        f_opt = self.optProb._mapObjtoOpt(f_user)
        assert_allclose(f_opt, f_user * self.objScale)
        assert_allclose(self.optProb._mapObjtoUser(f_opt), f_user)

    def test_mapCon_value_roundtrip(self):
        c_user = np.array([1.0, -2.0, 3.0])
        c_opt = self.optProb._mapContoOpt(c_user)
        assert_allclose(c_opt, c_user * self.conScale)
        assert_allclose(self.optProb._mapContoUser(c_opt), c_user)


class TestScalingEdgeCases(unittest.TestCase):
    def test_combined_scale_and_offset(self):
        """A DV group with both a non-unit scale and a non-zero offset is the
        classic place to get the order of operations wrong."""

        def objfunc(xdict):
            return {"obj": 0.0}, False

        optProb = Optimization("edge", objfunc)
        optProb.addVarGroup("x", 2, lower=-10, upper=10, scale=4.0, offset=3.0)
        optProb.addObj("obj")
        optProb.finalize()

        x_user = np.array([3.0, 7.0])  # note x_user[0] == offset
        x_opt = optProb._mapXtoOpt(x_user)
        # (x - 3) * 4
        assert_allclose(x_opt, np.array([0.0, 16.0]))
        assert_allclose(optProb._mapXtoUser(x_opt), x_user)

    def test_infinite_bounds_not_scaled(self):
        """INFINITY bounds must remain unbounded; scale/offset must not turn
        them into finite numbers in the assembled bounds."""
        # External modules
        from pyoptsparse.pyOpt_utils import INFINITY

        def objfunc(xdict):
            return {"obj": 0.0}, False

        optProb = Optimization("inf", objfunc)
        optProb.addVarGroup("x", 1, lower=None, upper=None, scale=10.0, offset=5.0)
        optProb.addObj("obj")
        optProb.finalize()

        var = optProb.variables["x"][0]
        # Variable stores scaled bounds; unbounded sides stay at exactly +/- INFINITY.
        self.assertEqual(var.lower, -INFINITY)
        self.assertEqual(var.upper, INFINITY)


if __name__ == "__main__":
    unittest.main()
