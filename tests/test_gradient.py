"""
Unit tests for the automatic-sensitivity engine (pyOpt_gradient.Gradient).

Previously only FD (the default) and CS were exercised, and only indirectly
through full optimizations where a derivative error shows up as slow/failed
convergence rather than a clear assertion. Here we drive the Gradient object
directly at a fixed point against a problem with a known analytic Jacobian, for
every sensitivity mode.

The test problem (all derivatives exact):
    obj = x0^2 + 2*x1^2 + 3*y^2
    c   = x0*x1 + y
with DV groups "x" (2 vars) and "y" (1 var).
"""

# Standard Python modules
import unittest

# External modules
import numpy as np
from numpy.testing import assert_allclose
from parameterized import parameterized

# First party modules
from pyoptsparse import Optimization
from pyoptsparse.pyOpt_gradient import Gradient

# Base point at which we evaluate the derivatives
X0 = {"x": np.array([1.5, -2.0]), "y": np.array([0.5])}

# Analytic Jacobian at X0 (user/unscaled space)
ANALYTIC = {
    "obj": {"x": np.array([3.0, -8.0]), "y": np.array([3.0])},
    "c": {"x": np.array([[-2.0, 1.5]]), "y": np.array([[1.0]])},
}

# Per-mode tolerances. Quadratic/bilinear functions are exact under central
# differencing and complex step; forward differencing carries an O(step) error.
TOLS = {
    "fd": dict(rtol=1e-4, atol=1e-5),
    "fdr": dict(rtol=1e-4, atol=1e-5),
    "cd": dict(rtol=1e-6, atol=1e-7),
    "cdr": dict(rtol=1e-6, atol=1e-7),
    "cs": dict(rtol=1e-11, atol=1e-12),
}


def objfunc(xdict):
    x = xdict["x"]
    y = xdict["y"]
    funcs = {}
    funcs["obj"] = x[0] ** 2 + 2 * x[1] ** 2 + 3 * y[0] ** 2
    funcs["c"] = np.array([x[0] * x[1] + y[0]])
    return funcs, False


def build_optProb(xScale=1.0, conScale=1.0):
    optProb = Optimization("grad-test", objfunc)
    optProb.addVarGroup("x", 2, lower=-10, upper=10, value=X0["x"], scale=xScale)
    optProb.addVarGroup("y", 1, lower=-10, upper=10, value=X0["y"], scale=xScale)
    optProb.addObj("obj")
    optProb.addConGroup("c", 1, lower=-100, upper=100, scale=conScale)
    optProb.finalize()
    return optProb


def assert_sens_matches_analytic(funcsSens, tol):
    for funcKey, perGroup in ANALYTIC.items():
        for dvGroup, expected in perGroup.items():
            assert_allclose(funcsSens[funcKey][dvGroup], expected, **tol)


class TestGradientModes(unittest.TestCase):
    @parameterized.expand(["fd", "fdr", "cd", "cdr", "cs"])
    def test_mode_matches_analytic(self, sensType):
        optProb = build_optProb()
        funcs, _ = objfunc(X0)
        grad = Gradient(optProb, sensType=sensType)
        funcsSens, fail = grad(X0, funcs)
        self.assertFalse(fail)
        assert_sens_matches_analytic(funcsSens, TOLS[sensType])

    def test_cs_returns_real(self):
        optProb = build_optProb()
        funcs, _ = objfunc(X0)
        grad = Gradient(optProb, sensType="cs")
        funcsSens, _ = grad(X0, funcs)
        for funcKey in ANALYTIC:
            for dvGroup in ANALYTIC[funcKey]:
                self.assertFalse(np.iscomplexobj(funcsSens[funcKey][dvGroup]))

    def test_fd_and_cs_agree(self):
        optProb = build_optProb()
        funcs, _ = objfunc(X0)
        fd = Gradient(optProb, sensType="fd")(X0, funcs)[0]
        cs = Gradient(optProb, sensType="cs")(X0, funcs)[0]
        for funcKey in ANALYTIC:
            for dvGroup in ANALYTIC[funcKey]:
                assert_allclose(fd[funcKey][dvGroup], cs[funcKey][dvGroup], rtol=1e-4, atol=1e-5)

    def test_default_step_sizes(self):
        # The defaults differ by mode; pin them so a refactor cannot silently
        # change differencing accuracy.
        self.assertEqual(Gradient(build_optProb(), "fd").sensStep, 1e-6)
        self.assertEqual(Gradient(build_optProb(), "fdr").sensStep, 1e-6)
        self.assertEqual(Gradient(build_optProb(), "cd").sensStep, 1e-4)
        self.assertEqual(Gradient(build_optProb(), "cdr").sensStep, 1e-4)
        self.assertEqual(Gradient(build_optProb(), "cs").sensStep, 1e-40j)


class TestRelativeStepAtZero(unittest.TestCase):
    """FDR/CDR use max(|step*x[i]|, step); at x[i]=0 this must fall back to the
    absolute floor rather than producing a zero step (which would divide by 0)."""

    @parameterized.expand(["fdr", "cdr"])
    def test_zero_dv(self, sensType):
        x0 = {"x": np.array([0.0, -2.0]), "y": np.array([0.5])}
        optProb = build_optProb()
        funcs, _ = objfunc(x0)
        grad = Gradient(optProb, sensType=sensType)
        funcsSens, fail = grad(x0, funcs)
        self.assertFalse(fail)
        # At x0=0: dobj/dx0 = 0, dc/dx0 = x1 = -2, dc/dx1 = x0 = 0
        self.assertTrue(np.all(np.isfinite(funcsSens["obj"]["x"])))
        assert_allclose(funcsSens["obj"]["x"], np.array([0.0, -8.0]), atol=1e-5)
        assert_allclose(funcsSens["c"]["x"], np.array([[-2.0, 0.0]]), atol=1e-5)


class TestGradientIgnoresScaling(unittest.TestCase):
    """CLAUDE.md: scaling must NOT be applied inside Gradient -- it operates on
    unscaled values and the scaling is reapplied later via the process* path.
    So adding DV/constraint scaling must not change the (unscaled) sens."""

    def test_scaling_does_not_affect_sens(self):
        optProb = build_optProb(xScale=7.0, conScale=0.3)
        funcs, _ = objfunc(X0)
        grad = Gradient(optProb, sensType="cs")
        funcsSens, _ = grad(X0, funcs)
        assert_sens_matches_analytic(funcsSens, TOLS["cs"])


if __name__ == "__main__":
    unittest.main()
