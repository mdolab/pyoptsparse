"""
Test uses Schittkowski's TP109 constraint problem.

    min 	3.0*x1+1.*10**(-6)*x1**3+0.522074*10**(-6)*x2**3+2.0*x2

    s.t.    -(x4-x3+0.550) <= 0
            -(x3-x4+0.550) <= 0
            -(2.25*10**(+6)-x1**2-x8**2) <= 0
            -(2.25*10**(+6)-x2**2-x9**2) <= 0
            -(x5*x6*np.sin(-x3-0.250) + +x5*x7*np.sin(-x4-0.250)+2.0*x5**2*b)*ra+400.0-x1 = 0
            -((x5*x6*np.sin(x3-0.250)+x6*x7*np.sin(x3-x4-0.250)+2.0*x6**2*b)*ra+400.0-x2) = 0
            -((x5*x7*np.sin(x4-0.250)+x6*x7*np.sin(x4-x3-0.250)+2.0*x7**2*b)*ra+881.7790) = 0
            -(x8+(x5*x6*np.cos(-x3-0.250)+x5*x7*np.cos(-x4-0.250)-2.0*x5**2*c)*ra+0.7533*10**(-3)*x5**2-200.00) = 0
            -(x9+(x5*x6*np.cos(x3-0.250)+x7*x6*np.cos(x3-x4-0.250)-2.0*x6**2*c)*ra+0.7533*10**(-3)*x(6)**2-200.00) = 0
            -((x5*x7*np.cos(x4-0.250)+x6*x7*np.cos(x4-x3-0.250)-2.0*x7**2*c)*ra+0.7533*10**(-3)*x7**2-22.9380) = 0
            0 <= xi, i = 1,2
            -0.55 <= xi <= 0.55, i = 3,4
            196.0 <= xi <= 252.0, i = 5,6,7
            -400.0 <= xi <= 800.0, i = 8,9

    where   a = 50.176
            b = np.sin(0.25)
            c = np.cos(0.25)
            ra = 1.0/50.176


    f*1 = 0.536206927538e+04
    x*1 = [0.674888100445e+03, 0.113417039470e+04, 0.133569060261e+00, -0.371152592466e+00, 0.252e+03, 0.252e+03, 0.201464535316e+03, 0.426660777226e+03, 0.368494083867e+03]
"""
# Standard Python modules
import unittest

# External modules
import numpy as np

# First party modules
from pyoptsparse import History, Optimization

# Local modules
from testing_utils import OptTest

USE_LINEAR = True


class TestTP109(OptTest):
    name = "TP109"
    DVs = {"xvars"}
    cons = {"con"}
    objs = {"obj"}
    fStar = 0.536206927538e4
    xStar = {
        "xvars": (
            675.0253392597418,
            1134.021089714841,
            0.1334850509182756,
            -0.37119026264778,
            252,
            252,
            201.4658577338128,
            426.6190056621182,
            368.4881990875219,
        )
    }

    def objfunc(self, xdict):
        x = xdict["xvars"]

        a = 50.1760
        b = np.sin(0.250)
        c = np.cos(0.250)
        funcs = {}
        funcs["obj"] = 3.0 * x[0] + (1e-6) * x[0] ** 3 + 0.522074e-6 * x[1] ** 3 + 2 * x[1]
        fcon = np.zeros(10, "D")
        fcon[0] = 2250000 - x[0] ** 2 - x[7] ** 2
        fcon[1] = 2250000 - x[1] ** 2 - x[8] ** 2
        fcon[2] = (
            x[4] * x[5] * np.sin(-x[2] - 0.25)
            + x[4] * x[6] * np.sin(-x[3] - 0.25)
            + 2 * b * x[4] ** 2
            - a * x[0]
            + 400 * a
        )
        fcon[3] = (
            x[4] * x[5] * np.sin(x[2] - 0.25)
            + x[5] * x[6] * np.sin(x[2] - x[3] - 0.25)
            + 2 * b * x[5] ** 2
            - a * x[1]
            + 400 * a
        )
        fcon[4] = (
            x[4] * x[6] * np.sin(x[3] - 0.25)
            + x[5] * x[6] * np.sin(x[3] - x[2] - 0.25)
            + 2 * b * x[6] ** 2
            + 881.779 * a
        )
        fcon[5] = (
            a * x[7]
            + x[4] * x[5] * np.cos(-x[2] - 0.25)
            + x[4] * x[6] * np.cos(-x[3] - 0.25)
            - 200 * a
            - 2 * c * x[4] ** 2
            + 0.7533e-3 * a * x[4] ** 2
        )
        fcon[6] = (
            a * x[8]
            + x[4] * x[5] * np.cos(x[2] - 0.25)
            + x[5] * x[6] * np.cos(x[2] - x[3] - 0.25)
            - 2 * c * x[5] ** 2
            + 0.7533e-3 * a * x[5] ** 2
            - 200 * a
        )
        fcon[7] = (
            x[4] * x[6] * np.cos(x[3] - 0.25)
            + x[5] * x[6] * np.cos(x[3] - x[2] - 0.25)
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

    def setup_optProb(self):
        # Optimization Object
        self.optProb = Optimization("TP109 Constraint Problem", self.objfunc)

        # Design Variables
        lower = [0.0, 0.0, -0.55, -0.55, 196, 196, 196, -400, -400]
        upper = [None, None, 0.55, 0.55, 252, 252, 252, 800, 800]
        value = [0, 0, 0, 0, 0, 0, 0, 0, 0]
        self.optProb.addVarGroup("xvars", 9, lower=lower, upper=upper, value=value)

        # Constraints
        lower = [0, 0, 0, 0, 0, 0, 0, 0]
        upper = [None, None, 0, 0, 0, 0, 0, 0]
        if not USE_LINEAR:
            lower.extend([0, 0])
            upper.extend([None, None])

        self.optProb.addConGroup("con", len(lower), lower=lower, upper=upper)

        # And the 2 linear constriants
        if USE_LINEAR:
            jac = np.zeros((1, 9))
            jac[0, 3] = 1.0
            jac[0, 2] = -1.0
            self.optProb.addConGroup(
                "lin_con", 1, lower=-0.55, upper=0.55, wrt=["xvars"], jac={"xvars": jac}, linear=True
            )

        # Objective
        self.optProb.addObj("obj")

    def test_snopt(self):
        self.optName = "SNOPT"
        self.setup_optProb()
        sol = self.optimize(storeHistory=True, sens="CS")
        self.assert_solution_allclose(sol, 1e-7)
        hist = History(self.histFileName)
        self.assertNotIn("lin_con", hist.getConNames())
        self.assertNotIn("lin_con", hist.getConInfo())
        val = hist.getValues()
        hist.getValues(scale=True)

        # Check that the function values in the history are real
        self.assertTrue(np.isrealobj(val["obj"]))
        self.assertTrue(np.isrealobj(val["con"]))

    def test_snopt_informs(self):
        self.optName = "SNOPT"
        self.setup_optProb()
        sol = self.optimize(optOptions={"Time Limit": 1e-15})
        self.assert_inform_equal(sol, 34)

    def test_slsqp(self):
        self.optName = "SLSQP"
        self.setup_optProb()
        sol = self.optimize(sens="CS")
        self.assert_solution_allclose(sol, 1e-7)
        # Check that the function values in the solution are real
        self.assertTrue(np.isrealobj(sol.objectives["obj"].value))
        self.assertTrue(np.isrealobj(sol.constraints["con"].value))

    def test_autorefine(self):
        self.optName = "ALPSO"
        self.setup_optProb()
        # Get first Solution
        sol1 = self.optimize(sens="CS")

        # Now re-optimize with SNOPT
        self.optProb = sol1
        self.optName = "SNOPT"
        self.setup_optProb()
        sol2 = self.optimize(sens="CS")

        # Check Solution
        self.assert_solution_allclose(sol2, 1e-2)


if __name__ == "__main__":
    unittest.main()
