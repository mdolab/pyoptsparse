"""Test class for Egor specific tests"""

import numpy as np

# First party modules
from pyoptsparse import OPT, Optimization
from pyoptsparse.testing import OptTest

import egobox as egx


class TestEgor(OptTest):

    def setup_xsinx_optProb(self):
        """
        Setup the optimization problem for the xsinx function.
        
        The xsinx function is defined as:
            f(x) = (x - 3.5) * sin((x - 3.5) / π)

        Domain: x ∈ [0, 25]
        Solution opt pb: fStar ≈ -15.1 at x ≈ 18.935
        """
        def objfunc(xdict):
            x = xdict["x"]
            funcs = {}
            # xsinx function: (x - 3.5) * sin((x - 3.5) / π)
            funcs["obj"] = (x - 3.5) * np.sin((x - 3.5) / np.pi)
            fail = False
            return funcs, fail

        optProb = Optimization("xsinx Function", objfunc)
        optProb.addVar("x", lower=0.0, upper=25.0)
        optProb.addObj("obj")
        self.optName = "Egor"
        self.optProb = optProb

    def test_egor(self):
        self.setup_xsinx_optProb()
        sol = self.optimize()
        # Check Solution
        print("xsinx Solution: ", sol.xStar, "f: ", sol.fStar)
        self.assertLess(sol.fStar, -15.1)  # Should find a negative minimum

    def test_egor_inform(self):
        self.setup_xsinx_optProb()
        # Test that the inform is "Maximum number of iterations reached"
        sol = self.optimize(optOptions={"max_iters": 1})
        self.assert_inform_equal(sol, int(egx.ExitStatus.MAX_ITERS_REACHED)) 

        # Test that the inform is "Target function value reached"
        sol = self.optimize(optOptions={"target": -10.})
        self.assert_inform_equal(sol, int(egx.ExitStatus.TARGET_COST_REACHED))  

        # Test that the inform is "Time limit reached"
        sol = self.optimize(optOptions={"timeout": 1e-6})
        self.assert_inform_equal(sol, int(egx.ExitStatus.TIMEOUT)) 

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
        print("Ackley Solution: ", sol.xStar, "f: ", sol.fStar)
        self.assertAlmostEqual(sol.fStar, 0.0, delta=1e-2)
        self.assertAlmostEqual(sol.xStar["xvars"][0], 0.0, delta=1e-2)
        self.assertAlmostEqual(sol.xStar["xvars"][1], 0.0, delta=1e-2)

