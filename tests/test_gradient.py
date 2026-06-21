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
X0 = {"x": [1.5, -2.0], "y": 0.5}

# Analytic Jacobian at X0
ANALYTIC = {
    "obj": {"x": [3.0, -8.0], "y": 3.0},
    "c": {"x": [[-2.0, 1.5]], "y": 1.0},
}


def objfunc(xdict):
    """
    Obj = x0^2 + 2*x1^2 + 3*y^2
    c   = x0*x1 + y
    """
    x, y = xdict["x"], xdict["y"]
    funcs = {}
    funcs["obj"] = x[0] ** 2 + 2 * x[1] ** 2 + 3 * y**2
    funcs["c"] = x[0] * x[1] + y
    return funcs, False


def build_optProb(objfun=objfunc, xScale=1.0, conScale=1.0):
    optProb = Optimization("grad-test", objfun)
    optProb.addVarGroup("x", 2, lower=-10, upper=10, value=X0["x"], scale=xScale)
    optProb.addVar("y", lower=-10, upper=10, value=X0["y"], scale=xScale)
    optProb.addObj("obj")
    optProb.addCon("c", lower=-100, upper=100, scale=conScale)
    optProb.finalize()
    return optProb


def assert_sens_matches_analytic(funcsSens, atol):
    for funcKey, perGroup in ANALYTIC.items():
        for dvGroup, expected in perGroup.items():
            assert_allclose(funcsSens[funcKey][dvGroup], expected, atol=atol)


class TestGradient(unittest.TestCase):
    @parameterized.expand(["fd", "fdr", "cd", "cdr", "cs"])
    def test_mode_matches_analytic(self, sensType):
        """Test that all differentiation modes produce Jacobians matching the analytic values,
        and that CS returns real (not complex) values."""
        optProb = build_optProb()
        funcs, _ = objfunc(X0)
        grad = Gradient(optProb, sensType=sensType)
        funcsSens, fail = grad(X0, funcs)
        self.assertFalse(fail)
        atol = 1e-12 if sensType == "cs" else 1e-5
        assert_sens_matches_analytic(funcsSens, atol=atol)

        if sensType == "cs":
            for funcKey in ANALYTIC:
                for dvGroup in ANALYTIC[funcKey]:
                    self.assertFalse(np.iscomplexobj(funcsSens[funcKey][dvGroup]))

    def test_failed_eval(self):
        """Test that failed flags from objfunc are caught by the Gradient class"""

        def always_fail(xdict):
            funcs, _ = objfunc(xdict)
            return funcs, True

        optProb = build_optProb(objfun=always_fail)
        funcs, _ = objfunc(X0)
        grad = Gradient(optProb, sensType="fd")
        _, fail = grad(X0, funcs)
        self.assertTrue(fail)

    def test_scaling(self):
        """Test that CS gradients are invariant under scaling"""
        optProb = build_optProb(xScale=7.0, conScale=0.3)
        funcs, _ = objfunc(X0)
        grad = Gradient(optProb, sensType="cs")
        funcsSens, _ = grad(X0, funcs)
        assert_sens_matches_analytic(funcsSens, 1e-12)


if __name__ == "__main__":
    unittest.main()
