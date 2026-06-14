"""Test class for Egor specific tests"""

import numpy as np

# First party modules
from pyoptsparse import OPT, Optimization
from pyoptsparse.testing import OptTest


class TestEgor(OptTest):

    def test_egor_ackley(self):
        """
        Test that Egor can optimize the Ackley function.
        """

        def objfunc(xdict):
            x = xdict["xvars"]
            funcs = {}
            funcs["obj"] = -20.0 * np.exp(-0.2 * np.sqrt(0.5 * (x[0] ** 2 + x[1] ** 2))) - np.exp(
                0.5 * (np.cos(2.0 * np.pi * x[0]) + np.cos(2.0 * np.pi * x[1]))
            ) + np.e + 20
            fail = False
            return funcs, fail

        optProb = Optimization("Ackley Function", objfunc)
        optProb.addVarGroup("xvars", 2, lower=[-32.768, -32.768], upper=[32.768, 32.768])
        optProb.addObj("obj")
        self.optName = "Egor"
        self.optProb = optProb
        sol = self.optimize(optOptions=
            {"max_iters": 100, 
             "verbose": 2, # level of verbosity, 0 = error, 1 = warn, 2 = info, 3 = debug
             "n_doe": 15, 
             "gp_config":{"corr_spec": 8}, # corr spec: 1 = absolute exponential, 2 = squared exponential, 4 = matern 3/2, 8 = matern 5/2
             "seed": 0})
        # Check Solution
        print("Solution: ", sol.xStar, "f: ", sol.fStar)
        self.assertAlmostEqual(sol.fStar, 0.0, delta=1e-2)
        self.assertAlmostEqual(sol.xStar["xvars"][0], 0.0, delta=1e-2)
        self.assertAlmostEqual(sol.xStar["xvars"][1], 0.0, delta=1e-2)
