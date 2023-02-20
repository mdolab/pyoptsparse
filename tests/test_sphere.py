"""Test solution of Sphere problem"""

# Standard Python modules
import unittest

# External modules
import numpy as np
from parameterized import parameterized

# First party modules
from pyoptsparse import Optimization

# Local modules
from testing_utils import OptTest


class TestSphere(OptTest):
    ## Solve unconstrained Sphere problem.
    #  This problem is scalable w.r.t. design variables number.
    #  We select a problem with 4 design variables, but the
    #  location and value of the minimum do not change with DV
    #  dimensionality
    #
    #
    #  min   Sum(x[i] ** 2)
    #
    #  The minimum is located at x=(0,....,0) where x
    #  is an arbitrarily sized vector depending on the number N
    #  of design variables.
    #  At the optimum, the function is f(x) = 0.
    #  We select a random initial point for our test.
    ##

    name = "Sphere"
    N = 4
    objs = {"obj"}
    cons = set()
    DVs = {"xvars"}
    fStar = 0.0
    xStar = {"xvars": np.zeros(N)}

    # Tolerances
    tol = {"ALPSO": 1e-3}

    optOptions = {
        "ALPSO": {  # sphere
            "SwarmSize": 20,
            "maxOuterIter": 10,
            "c1": 1.0,  # Cognitive Parameter
            "c2": 1.25,  # Social Parameter
            "stopCriteria": 0,  # 0: maxOuterIter, 1: convergence
            "seed": 1235,
        }
    }

    def objfunc(self, xdict):
        self.nf += 1
        x = xdict["xvars"]

        funcs = {"obj": np.dot(x, x)}

        fail = False
        return funcs, fail

    def sens(self, xdict, _funcs):
        self.ng += 1
        x = xdict["xvars"]

        funcsSens = {"obj": {"xvars": 2 * np.ones(len(x))}}

        fail = False
        return funcsSens, fail

    def setup_optProb(self):
        # Optimization Object
        self.optProb = Optimization("Sphere Problem", self.objfunc)

        np.random.seed(10)
        value = np.random.normal(size=self.N)

        lower = np.ones(self.N) * -50
        upper = np.ones(self.N) * 50
        self.optProb.addVarGroup("xvars", self.N, lower=lower, upper=upper, value=value)

        # Objective
        self.optProb.addObj("obj")

    @parameterized.expand(["ALPSO"])
    def test_optimization(self, optName):
        self.optName = optName
        self.setup_optProb()
        optOptions = self.optOptions.pop(optName, None)
        self.optimize_with_hotstart(self.tol[optName], optOptions=optOptions)


if __name__ == "__main__":
    unittest.main()
