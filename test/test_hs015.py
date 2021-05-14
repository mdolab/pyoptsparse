"""Test solution of problem HS15 from the Hock & Schittkowski collection"""

# Standard Python modules
import unittest

# External modules
import numpy as np
from parameterized import parameterized

# First party modules
from pyoptsparse import History, Optimization

# Local modules
from testing_utils import OptTest


class TestHS15(OptTest):

    ## Solve test problem HS15 from the Hock & Schittkowski collection.
    #
    #  min   100 (x2 - x1^2)^2 + (1 - x1)^2
    #  s.t.  x1 x2 >= 1
    #        x1 + x2^2 >= 0
    #        x1 <= 0.5
    #
    #  The standard start point (-2, 1) usually converges to the standard
    #  minimum at (0.5, 2.0), with final objective = 306.5.
    #  Sometimes the solver converges to another local minimum
    #  at (-0.79212, -1.26243), with final objective = 360.4.
    ##

    name = "HS015"
    DVs = {"xvars"}
    cons = {"con"}
    objs = {"obj"}
    extras = {"extra1", "extra2"}
    fStar = [
        306.5,
        360.379767,
    ]
    xStar = [
        {"xvars": (0.5, 2.0)},
        {"xvars": (-0.79212322, -1.26242985)},
    ]
    tol = {
        "SLSQP": 1e-8,
        "NLPQLP": 1e-12,
        "IPOPT": 1e-4,
        "ParOpt": 1e-6,
        "CONMIN": 1e-10,
        "PSQP": 5e-12,
    }
    optOptions = {}

    def objfunc(self, xdict):
        self.nf += 1
        x = xdict["xvars"]
        funcs = {}
        funcs["obj"] = [100 * (x[1] - x[0] ** 2) ** 2 + (1 - x[0]) ** 2]
        conval = np.zeros(2, "D")
        conval[0] = x[0] * x[1]
        conval[1] = x[0] + x[1] ** 2
        funcs["con"] = conval
        # extra keys
        funcs["extra1"] = 0.0
        funcs["extra2"] = 1.0
        fail = False
        return funcs, fail

    def sens(self, xdict, funcs):
        self.ng += 1
        x = xdict["xvars"]
        funcsSens = {}
        funcsSens["obj"] = {
            "xvars": [2 * 100 * (x[1] - x[0] ** 2) * (-2 * x[0]) - 2 * (1 - x[0]), 2 * 100 * (x[1] - x[0] ** 2)]
        }
        funcsSens["con"] = {"xvars": [[x[1], x[0]], [1, 2 * x[1]]]}
        fail = False
        return funcsSens, fail

    def setup_optProb(self):
        # Optimization Object
        self.optProb = Optimization("HS15 Constraint Problem", self.objfunc)

        # Design Variables
        lower = [-5.0, -5.0]
        upper = [0.5, 5.0]
        value = [-2, 1.0]
        self.optProb.addVarGroup("xvars", 2, lower=lower, upper=upper, value=value)

        # Constraints
        lower = [1.0, 0.0]
        upper = [None, None]
        self.optProb.addConGroup("con", 2, lower=lower, upper=upper)

        # Objective
        self.optProb.addObj("obj")

    def test_snopt(self):
        self.optName = "SNOPT"
        self.setup_optProb()
        store_vars = ["step", "merit", "feasibility", "optimality", "penalty", "Hessian", "condZHZ", "slack", "lambda"]
        optOptions = {"Save major iteration variables": store_vars}
        self.optimize_with_hotstart(1e-12, optOptions=optOptions)

        hist = History(self.histFileName, flag="r")
        data = hist.getValues(callCounters=["last"])
        keys = hist.getIterKeys()
        self.assertIn("isMajor", keys)
        self.assertEqual(7, data["nMajor"])
        for var in store_vars:
            self.assertIn(var, data.keys())
        self.assertEqual(data["Hessian"].shape, (1, 2, 2))
        self.assertEqual(data["feasibility"].shape, (1, 1))
        self.assertEqual(data["slack"].shape, (1, 2))
        self.assertEqual(data["lambda"].shape, (1, 2))
        # dv = sol.getDVs()
        # sol_xvars = [sol.variables["xvars"][i].value for i in range(2)]
        # assert_allclose(sol_xvars, dv["xvars"], atol=tol, rtol=tol)

    @parameterized.expand(["IPOPT", "SLSQP", "PSQP", "CONMIN", "NLPQLP", "ParOpt"])
    def test_optimization(self, optName):
        self.optName = optName
        self.setup_optProb()
        optOptions = self.optOptions.pop(optName, None)
        sol = self.optimize(optOptions=optOptions)
        # Check Solution
        self.assert_solution_allclose(sol, self.tol[optName])
        # Check informs
        self.assert_inform_equal(sol)


if __name__ == "__main__":
    unittest.main()
