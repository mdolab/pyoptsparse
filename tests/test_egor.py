"""Test class for Egor specific tests"""

import numpy as np
import tempfile

# First party modules
from pyoptsparse import Optimization
from pyoptsparse.testing import OptTest


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
        self.assert_inform_equal(sol, 1)

        # Test that the inform is "Target function value reached"
        sol = self.optimize(optOptions={"target": -10.0})
        self.assert_inform_equal(sol, 2)

        # Test that the inform is "Time limit reached"
        sol = self.optimize(optOptions={"timeout": 1e-6})
        self.assert_inform_equal(sol, 5)

    def test_egor_warm_start(self):
        with tempfile.TemporaryDirectory() as outdir:
            self.setup_xsinx_optProb()
            # First run to generate a history file
            sol1 = self.optimize(optOptions={"max_iters": 1, "outdir": outdir, "seed": 0})
            print("First run: ", sol1.xStar, "f: ", sol1.fStar)
            # Second run with warm start
            sol2 = self.optimize(optOptions={"max_iters": 5, "outdir": outdir, "warm_start": True})
            print("Second run (warm start): ", sol2.xStar, "f: ", sol2.fStar)
            # Check that the second run continued from the first run
            self.assertGreater(sol1.fStar, sol2.fStar)

    def test_egor_config(self):
        with tempfile.TemporaryDirectory() as outdir:
            self.setup_xsinx_optProb()
            # Test that the gp_config option is passed correctly
            gp_config = {"corr_spec": 4, "kpls_dim": 1}
            _ = self.optimize(
                optOptions={
                    "infill_strategy": 1,
                    "gp_config": gp_config,
                    "outdir": outdir,
                    "trego": {"n_gl_steps": (1, 3)},
                }
            )
            # read egor_config.json from outdir and check that corr_spec is 4
            import json

            with open(f"{outdir}/egor_config.json", "r") as f:
                egor_config = json.load(f)
            print("Egor config: ", egor_config)
            self.assertEqual(egor_config["gp"]["correlation_spec"], "MATERN32")
            self.assertEqual(egor_config["gp"]["kpls_dim"], 1)
            self.assertEqual(egor_config["infill_criterion"]["type_infill"], "ExpectedImprovement")
            self.assertEqual(egor_config["iteration_strategy"]["type_iteration_strategy"], "TregoStrategy")
            self.assertEqual(egor_config["iteration_strategy"]["n_gl_steps"], [1, 3])

    def test_egor_ackley(self):
        """
        Test that Egor can optimize the Ackley function.
        """

        def objfunc(xdict):
            x = xdict["xvars"]
            funcs = {}
            funcs["obj"] = (
                -20.0 * np.exp(-0.2 * np.sqrt(0.5 * (x[0] ** 2 + x[1] ** 2)))
                - np.exp(0.5 * (np.cos(2.0 * np.pi * x[0]) + np.cos(2.0 * np.pi * x[1])))
                + np.e
                + 20
            )
            fail = False
            return funcs, fail

        optProb = Optimization("Ackley Function", objfunc)
        optProb.addVarGroup("xvars", 2, lower=[-32.768, -32.768], upper=[32.768, 32.768])
        optProb.addObj("obj")
        self.optName = "Egor"
        self.optProb = optProb
        sol = self.optimize(
            optOptions={
                "max_iters": 100,
                "verbose": 2,  # level of verbosity, 0 = error, 1 = warn, 2 = info, 3 = debug
                "n_doe": 15,
                "gp_config": {
                    "corr_spec": 8
                },  # corr spec: 1 = absolute exponential, 2 = squared exponential, 4 = matern 3/2, 8 = matern 5/2
                "seed": 0,
            }
        )
        # Check Solution
        print("Ackley Solution: ", sol.xStar, "f: ", sol.fStar)
        self.assertAlmostEqual(sol.fStar, 0.0, delta=1e-2)
        self.assertAlmostEqual(sol.xStar["xvars"][0], 0.0, delta=1e-2)
        self.assertAlmostEqual(sol.xStar["xvars"][1], 0.0, delta=1e-2)

    def test_egor_g24(self):
        """
        Test Egor on the G24 problem.

        The G24 problem is defined as:
            minimize f(x) = -x1 - x2
            subject to:
                c1(x) = -2*x1^4 + 8*x1^3 - 8*x1^2 + x2 - 2 <= 0
                c2(x) = -4*x1^4 + 32*x1^3 - 88*x1^2 + 96*x1 + x2 - 36 <= 0
            with x1 in [0, 3] and x2 in [0, 4]

        Global optimum: x_opt = (2.3295, 3.1785), f_opt = -5.5080
        """

        def objfunc(xdict):
            x = xdict["xvars"]
            funcs = {}
            funcs["obj"] = -x[0] - x[1]
            funcs["con"] = [
                -2.0 * x[0] ** 4 + 8.0 * x[0] ** 3 - 8.0 * x[0] ** 2 + x[1] - 2.0,
                -4.0 * x[0] ** 4 + 32.0 * x[0] ** 3 - 88.0 * x[0] ** 2 + 96.0 * x[0] + x[1] - 36.0,
            ]
            fail = False
            return funcs, fail

        optProb = Optimization("G24 Function", objfunc)
        optProb.addVarGroup("xvars", 2, lower=[0.0, 0.0], upper=[3.0, 4.0])
        optProb.addObj("obj")
        optProb.addConGroup("con", 2, upper=0.0)
        self.optName = "Egor"
        self.optProb = optProb
        sol = self.optimize(
            optOptions={"max_iters": 30, "n_doe": 5, "target": -5.50, "cstr_tol": [1e-3, 1e-3], "verbose": 2}
        )
        # Check Solution
        print("G24 Solution: ", sol.xStar, "f: ", sol.fStar)
        self.assertLess(sol.fStar, -5.50)  # Should find a value close to -5.5080
        self.assertAlmostEqual(sol.xStar["xvars"][0], 2.3295, delta=0.1)
        self.assertAlmostEqual(sol.xStar["xvars"][1], 3.1785, delta=0.1)
