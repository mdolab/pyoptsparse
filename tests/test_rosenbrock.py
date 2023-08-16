"""Test solution of Rosenbrock problem"""

# Standard Python modules
import sys
import unittest

# External modules
import numpy as np
from parameterized import parameterized
from sqlitedict import SqliteDict

# First party modules
from pyoptsparse import History, Optimization

# Local modules
from testing_utils import OptTest


class TestRosenbrock(OptTest):
    ## Solve unconstrained Rosenbrock problem.
    #  This problem is scalable w.r.t. design variables number.
    #  We select a problem with 4 design variables, but the
    #  location and value of the minimum do not change with DV
    #  dimensionality
    #
    #
    #  min   100 * (x[i + 1] - x[i] ** 2) ** 2 + (1 - x[i]) ** 2
    #
    #  The minimum is located at x=(1,....,1) where x
    #  is an arbitrarily sized vector depending on the number N
    #  of design variables.
    #  At the optimum, the function is f(x) = 0.
    #  We select a random initial point for our test.
    ##

    name = "Rosenbrock"
    N = 4
    objs = {"obj"}
    cons = set()
    DVs = {"xvars"}
    fStar = 0.0
    xStar = {"xvars": np.ones(N)}

    # Tolerances
    tol = {
        "SNOPT": 1e-6,
        "IPOPT": 1e-6,
        "NLPQLP": 1e-6,
        "SLSQP": 1e-6,
        "CONMIN": 1e-9,
        "PSQP": 1e-8,
        "ParOpt": 1e-8,
    }
    optOptions = {
        "SLSQP": {"ACC": 1e-10},
        "NLPQLP": {"accuracy": 1e-10},
    }

    def objfunc(self, xdict):
        self.nf += 1
        x = xdict["xvars"]

        funcs = {}
        funcs["obj"] = 0

        for i in range(len(x) - 1):
            funcs["obj"] += 100 * (x[i + 1] - x[i] ** 2) ** 2 + (1 - x[i]) ** 2

        fail = False
        return funcs, fail

    def sens(self, xdict, funcs):
        self.ng += 1
        x = xdict["xvars"]
        funcsSens = {}
        grads = np.zeros(len(x))

        for i in range(len(x) - 1):
            grads[i] += 2 * (200 * x[i] ** 3 - 200 * x[i] * x[i + 1] + x[i] - 1)
            grads[i + 1] += 200 * (x[i + 1] - x[i] ** 2)

        funcsSens["obj"] = {"xvars": grads}

        fail = False
        return funcsSens, fail

    def setup_optProb(self):
        # Optimization Object
        self.optProb = Optimization("Rosenbrock Problem", self.objfunc)

        np.random.seed(10)
        value = np.random.normal(size=self.N)

        lower = np.ones(self.N) * -50
        upper = np.ones(self.N) * 50
        self.optProb.addVarGroup("xvars", self.N, lower=lower, upper=upper, value=value)

        # Objective
        self.optProb.addObj("obj")

    def test_snopt(self):
        self.optName = "SNOPT"
        self.setup_optProb()
        store_vars = ["Hessian", "slack", "lambda", "penalty_vector", "nS", "BSwap", "maxVi"]
        optOptions = {"Save major iteration variables": store_vars}
        self.optimize_with_hotstart(1e-8, optOptions=optOptions)

        hist = History(self.histFileName, flag="r")
        data = hist.getValues(callCounters=["last"])
        keys = hist.getIterKeys()
        self.assertIn("isMajor", keys)
        self.assertEqual(36, data["nMajor"])
        for var in store_vars:
            self.assertIn(var, data.keys())
        self.assertEqual(data["Hessian"].shape, (1, 4, 4))
        self.assertEqual(data["feasibility"].shape, (1, 1))
        self.assertEqual(data["slack"].shape, (1, 1))
        self.assertEqual(data["lambda"].shape, (1, 1))

    def test_snopt_hotstart_starting_from_grad(self):
        self.optName = "SNOPT"
        self.setup_optProb()
        histName = f"{self.id()}.hst"

        # Optimize without hot start and store the history
        self.optimize(storeHistory=histName, hotStart=False)

        # Load the history dictionary
        hist = SqliteDict(histName)

        # Delete the last two keys in the dictionary
        lastKey = hist["last"]
        for i in range(2):
            del hist[str(int(lastKey) - i)]
        hist.commit()

        # Optimize starting from the modified history file
        # The first call will be a gradient evaluation
        self.optimize(storeHistory=False, hotStart=histName)

        # Check that we had two function evaluations
        # The first is from a recursive call and the second is the 'last' call we deleted
        self.assertEqual(self.nf, 2)

        # Also check that we had two gradient evaluations
        # The first is from a call we deleted and the second is the call after 'last'
        self.assertEqual(self.ng, 2)

    @parameterized.expand(["IPOPT", "SLSQP", "PSQP", "CONMIN", "NLPQLP", "ParOpt"])
    def test_optimization(self, optName):
        self.optName = optName
        if optName == "IPOPT" and sys.platform == "win32":
            raise unittest.SkipTest()
        self.setup_optProb()
        optOptions = self.optOptions.pop(optName, None)
        self.optimize_with_hotstart(self.tol[optName], optOptions=optOptions)


if __name__ == "__main__":
    unittest.main()
