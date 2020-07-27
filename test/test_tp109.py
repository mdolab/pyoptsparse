"""
Test uses Schittkowski's TP109 constraint problem.

    min 	3.0*x1+1.*10**(-6)*x1**3+0.522074*10**(-6)*x2**3+2.0*x2
    s.t.    -(x4-x3+0.550) <= 0
            -(x3-x4+0.550) <= 0
            -(2.25*10**(+6)-x1**2-x8**2) <= 0
            -(2.25*10**(+6)-x2**2-x9**2) <= 0
            -(x5*x6*sin(-x3-0.250) + +x5*x7*sin(-x4-0.250)+2.0*x5**2*b)*ra+400.0-x1 = 0
            -((x5*x6*sin(x3-0.250)+x6*x7*sin(x3-x4-0.250)+2.0*x6**2*b)*ra+400.0-x2) = 0
            -((x5*x7*sin(x4-0.250)+x6*x7*sin(x4-x3-0.250)+2.0*x7**2*b)*ra+881.7790) = 0
            -(x8+(x5*x6*cos(-x3-0.250)+x5*x7*cos(-x4-0.250)-2.0*x5**2*c)*ra+0.7533*10**(-3)*x5**2-200.00) = 0
            -(x9+(x5*x6*cos(x3-0.250)+x7*x6*cos(x3-x4-0.250)-2.0*x6**2*c)*ra+0.7533*10**(-3)*x(6)**2-200.00) = 0
            -((x5*x7*cos(x4-0.250)+x6*x7*cos(x4-x3-0.250)-2.0*x7**2*c)*ra+0.7533*10**(-3)*x7**2-22.9380) = 0
            0 <= xi, i = 1,2
            -0.55 <= xi <= 0.55, i = 3,4
            196.0 <= xi <= 252.0, i = 5,6,7
            -400.0 <= xi <= 800.0, i = 8,9
    where   a = 50.176
            b = sin(0.25)
            c = cos(0.25)
            ra = 1.0/50.176


    f*1 = 0.536206927538e+04
    x*1 = [0.674888100445e+03, 0.113417039470e+04, 0.133569060261e+00, -0.371152592466e+00, 0.252e+03, 0.252e+03, 0.201464535316e+03, 0.426660777226e+03, 0.368494083867e+03]
"""
import unittest

import numpy as np
from numpy import sin, cos
from numpy.testing import assert_allclose
from pyoptsparse import Optimization, OPT, History
from pyoptsparse.pyOpt_error import Error

USE_LINEAR = True


def objfunc(xdict):
    x = xdict["xvars"]

    a = 50.1760
    b = sin(0.250)
    c = cos(0.250)
    funcs = {}
    funcs["obj"] = 3.0 * x[0] + (1e-6) * x[0] ** 3 + 0.522074e-6 * x[1] ** 3 + 2 * x[1]
    fcon = np.zeros(10, "D")
    fcon[0] = 2250000 - x[0] ** 2 - x[7] ** 2
    fcon[1] = 2250000 - x[1] ** 2 - x[8] ** 2
    fcon[2] = x[4] * x[5] * sin(-x[2] - 0.25) + x[4] * x[6] * sin(-x[3] - 0.25) + 2 * b * x[4] ** 2 - a * x[0] + 400 * a
    fcon[3] = (
        x[4] * x[5] * sin(x[2] - 0.25) + x[5] * x[6] * sin(x[2] - x[3] - 0.25) + 2 * b * x[5] ** 2 - a * x[1] + 400 * a
    )
    fcon[4] = x[4] * x[6] * sin(x[3] - 0.25) + x[5] * x[6] * sin(x[3] - x[2] - 0.25) + 2 * b * x[6] ** 2 + 881.779 * a
    fcon[5] = (
        a * x[7]
        + x[4] * x[5] * cos(-x[2] - 0.25)
        + x[4] * x[6] * cos(-x[3] - 0.25)
        - 200 * a
        - 2 * c * x[4] ** 2
        + 0.7533e-3 * a * x[4] ** 2
    )
    fcon[6] = (
        a * x[8]
        + x[4] * x[5] * cos(x[2] - 0.25)
        + x[5] * x[6] * cos(x[2] - x[3] - 0.25)
        - 2 * c * x[5] ** 2
        + 0.7533e-3 * a * x[5] ** 2
        - 200 * a
    )
    fcon[7] = (
        x[4] * x[6] * cos(x[3] - 0.25)
        + x[5] * x[6] * cos(x[3] - x[2] - 0.25)
        - 2 * c * x[6] ** 2
        - 22.938 * a
        + 0.7533e-3 * a * x[6] ** 2
    )
    fcon[8] = x[3] - x[2] + 0.55
    fcon[9] = x[2] - x[3] + 0.55

    if USE_LINEAR:
        funcs["con"] = fcon[0:8]
    else:
        funcs["con"] = fcon[0:10]
    fail = False

    return funcs, fail


class TestTP109(unittest.TestCase):
    def optimize(self, optName, tol, optOptions={}, storeHistory=None):
        # Optimization Object
        optProb = Optimization("TP109 Constraint Problem", objfunc)

        # Design Variables
        lower = [0.0, 0.0, -0.55, -0.55, 196, 196, 196, -400, -400]
        upper = [None, None, 0.55, 0.55, 252, 252, 252, 800, 800]
        value = [0, 0, 0, 0, 0, 0, 0, 0, 0]
        optProb.addVarGroup("xvars", 9, lower=lower, upper=upper, value=value)

        # Constraints
        lower = [0, 0, 0, 0, 0, 0, 0, 0]
        upper = [None, None, 0, 0, 0, 0, 0, 0]
        if not USE_LINEAR:
            lower.extend([0, 0])
            upper.extend([None, None])

        optProb.addConGroup("con", len(lower), lower=lower, upper=upper)

        # And the 2 linear constriants
        if USE_LINEAR:
            jac = np.zeros((1, 9))
            jac[0, 3] = 1.0
            jac[0, 2] = -1.0
            optProb.addConGroup("lin_con", 1, lower=-0.55, upper=0.55, wrt=["xvars"], jac={"xvars": jac}, linear=True)

        # Objective
        optProb.addObj("obj")

        # Check optimization problem:
        # optProb.printSparsity()

        # Optimizer
        try:
            opt = OPT(optName, options=optOptions)
        except Error:
            raise unittest.SkipTest("Optimizer not available:", optName)

        # Solution
        sol = opt(optProb, sens="CS", storeHistory=storeHistory)

        # Check Solution
        assert_allclose(sol.objectives["obj"].value, 0.536206927538e04, atol=tol, rtol=tol)

    def test_snopt(self):
        name = "tp109_snopt.hst"
        self.optimize("snopt", 1e-7, storeHistory=name)
        hist = History(name)
        self.assertNotIn("lin_con", hist.getConNames())
        self.assertNotIn("lin_con", hist.getConInfo())
        hist.getValues()
        hist.getValues(scale=True)

    def test_slsqp(self):
        self.optimize("slsqp", 1e-7)

    def test_autorefine(self):
        # Optimization Object
        optProb = Optimization("TP109 Constraint Problem", objfunc)

        # Design Variables (Removed infinite bounds for ALPSO)
        lower = [0.0, 0.0, -0.55, -0.55, 196, 196, 196, -400, -400]
        upper = [2000, 2000, 0.55, 0.55, 252, 252, 252, 800, 800]
        value = [0, 0, 0, 0, 0, 0, 0, 0, 0]
        optProb.addVarGroup("xvars", 9, lower=lower, upper=upper, value=value)

        # Constraints
        lower = [0, 0, 0, 0, 0, 0, 0, 0]
        upper = [None, None, 0, 0, 0, 0, 0, 0]
        if not USE_LINEAR:
            lower.extend([0, 0])
            upper.extend([None, None])

        optProb.addConGroup("con", len(lower), lower=lower, upper=upper)

        # And the 2 linear constriants
        if USE_LINEAR:
            jac = np.zeros((1, 9))
            jac[0, 3] = 1.0
            jac[0, 2] = -1.0
            optProb.addConGroup("lin_con", 1, lower=-0.55, upper=0.55, wrt=["xvars"], jac={"xvars": jac}, linear=True)

        # Objective
        optProb.addObj("obj")

        # Check optimization problem:
        # optProb.printSparsity()

        # Global Optimizer: ALPSO
        try:
            opt1 = OPT("ALPSO")
        except Error:
            raise unittest.SkipTest("Optimizer not available:", "ALPSO")

        # Get first Solution
        sol1 = opt1(optProb)

        # Now run the previous solution with SNOPT
        try:
            opt2 = OPT("SNOPT")
        except Error:
            raise unittest.SkipTest("Optimizer not available:", "SNOPT")

        sol2 = opt2(sol1)

        # Check Solution
        assert_allclose(sol2.objectives["obj"].value, 0.536206927538e04, atol=1e-2, rtol=1e-2)


if __name__ == "__main__":
    unittest.main()
