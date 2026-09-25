"""Test solution of problem HS60 from the Hock & Schittkowski collection"""

# Standard Python modules
import unittest

# External modules
import numpy as np
from parameterized import parameterized

# First party modules
from pyoptsparse import Optimization
from pyoptsparse.testing import OptTest


class TestHS60(OptTest):
    ## Solve test problem HS60 from the Hock & Schittkowski collection.
    #
    #  min   (x1 - 1)^2 + (x1 - x2)^2 + (x2 - x3)^4
    #  s.t.  x1 (1 + x2^2) + x3^4 = 4 + 3 sqrt(2)
    #        -10 <= xi <= 10

    name = "HS060"
    DVs = {"xvars"}
    cons = {"con"}
    objs = {"obj"}
    fStar = 0.0325682002513
    xStar = {"xvars": (1.10485902423, 1.19667419413, 1.53526225739)}

    tol = {
        "SNOPT": 1e-6,
        "IPOPT": 1e-6,
        "SLSQP": 1e-6,
        "PSQP": 1e-6,
        "NLPQLP": 3e-4,
        "Uno": 1e-4,
        "Egor": 7e-2,
    }
    optOptions = {
        "Egor": {"max_iters": 100, "n_doe": 20, "seed": 42},
    }

    def objfunc(self, xdict):
        x = xdict["xvars"]
        funcs = {}
        funcs["obj"] = (x[0] - 1) ** 2 + (x[0] - x[1]) ** 2 + (x[1] - x[2]) ** 4
        funcs["con"] = [x[0] * (1 + x[1] ** 2) + x[2] ** 4]
        fail = False
        return funcs, fail

    def sens(self, xdict, funcs):
        x = xdict["xvars"]
        funcsSens = {}
        funcsSens["obj"] = {
            "xvars": [
                2 * (x[0] - 1) + 2 * (x[0] - x[1]),
                -2 * (x[0] - x[1]) + 4 * (x[1] - x[2]) ** 3,
                -4 * (x[1] - x[2]) ** 3,
            ]
        }
        funcsSens["con"] = {"xvars": [[1 + x[1] ** 2, 2 * x[0] * x[1], 4 * x[2] ** 3]]}
        fail = False
        return funcsSens, fail

    def setup_optProb(self):
        self.optProb = Optimization("HS60 Constraint Problem", self.objfunc, sens=self.sens)
        self.optProb.addVarGroup("xvars", 3, lower=-10, upper=10, value=2.0)
        rhs = 4 + 3 * np.sqrt(2)
        self.optProb.addConGroup("con", 1, lower=rhs, upper=rhs)
        self.optProb.addObj("obj")

    @parameterized.expand(["SNOPT", "IPOPT", "SLSQP", "PSQP", "NLPQLP", "Uno", "Egor"])
    def test_optimization(self, optName):
        self.optName = optName
        self.setup_optProb()
        optOptions = self.optOptions.get(optName, None)
        sol = self.optimize(optOptions=optOptions)
        self.assert_solution_allclose(sol, self.tol[optName])
        self.assert_inform_equal(sol)


if __name__ == "__main__":
    unittest.main()
