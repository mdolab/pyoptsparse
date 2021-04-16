"""Test solution of problem HS71 from the Hock & Schittkowski collection"""

# Standard Python modules
import unittest

# External modules
import numpy as np
from numpy.testing import assert_allclose

# First party modules
from pyoptsparse import History, Optimization
from utils import OptTest
from parameterized import parameterized


class TestHS71(OptTest):
    # Optimization problem definition
    name = "hs071"
    fStar = 17.0140172
    xStar = [1.0, 4.743, 3.82115, 1.37941]
    lambdaStar = [0.55229366, -0.16146857]

    # Optimizer
    optOptions = {}
    tol = {
        "SNOPT": 1e-6,
        "IPOPT": 1e-6,
        "NLPQLP": 1e-6,
        "SLSQP": 1e-6,
        "CONMIN": 1e-2,
        "PSQP": 1e-6,
        "ParOpt": 1e-6,
    }

    def objfunc(self, xdict):
        x = xdict["xvars"]
        funcs = {}
        funcs["obj"] = x[0] * x[3] * (x[0] + x[1] + x[2]) + x[2]
        funcs["con"] = [x[0] * x[1] * x[2] * x[3], x[0] * x[0] + x[1] * x[1] + x[2] * x[2] + x[3] * x[3]]
        fail = False
        return funcs, fail

    def sens(self, xdict, funcs):
        x = xdict["xvars"]
        funcsSens = {}
        funcsSens["obj"] = {
            "xvars": np.array(
                [x[0] * x[3] + x[3] * (x[0] + x[1] + x[2]), x[0] * x[3], x[0] * x[3] + 1.0, x[0] * (x[0] + x[1] + x[2])]
            )
        }
        jac = [
            [x[1] * x[2] * x[3], x[0] * x[2] * x[3], x[0] * x[1] * x[3], x[0] * x[1] * x[2]],
            [2.0 * x[0], 2.0 * x[1], 2.0 * x[2], 2.0 * x[3]],
        ]
        funcsSens["con"] = {"xvars": jac}
        fail = False
        return funcsSens, fail

    def setup_optProb(self, xScale=1.0, objScale=1.0, conScale=1.0, offset=0.0, check_solution=True):
        # Optimization Object
        optProb = Optimization("HS071 Constraint Problem", self.objfunc, sens=self.sens)

        # Design Variables
        x0 = [1.0, 5.0, 5.0, 1.0]
        optProb.addVarGroup("xvars", 4, lower=1, upper=5, value=x0, scale=xScale, offset=offset)

        # Constraints
        optProb.addConGroup("con", 2, lower=[25, 40], upper=[None, 40], scale=conScale)

        # Objective
        optProb.addObj("obj", scale=objScale)
        return optProb

    def test_slsqp_setDV(self):
        """
        Test that setDV works as expected, even with scaling/offset
        """
        test_name = "hs071_SLSQP_setDV"
        histFileName = "{}.hst".format(test_name)
        newDV = {"xvars": np.array([1, 4, 4, 1])}
        optOptions = {"IFILE": "{}.out".format(test_name)}
        optProb = self.setup_optProb(xScale=1.5, conScale=1.2, objScale=32, offset=1.5)
        sol = self.optimize(optProb, "SLSQP", setDV=newDV, storeHistory=histFileName, optOptions=optOptions)
        self.check_solution(sol, self.tol["SLSQP"])
        # Verify the history file
        hist = History(histFileName, flag="r")
        init = hist.getValues(names="xvars", callCounters="0", scale=False)
        x_init = init["xvars"][0]
        assert_allclose(x_init, newDV["xvars"], atol=1e-5, rtol=1e-5)

    def test_snopt_setDVFromHist(self):
        """
        Test that setDVFromHistory works as expected, even with scaling/offset
        """
        test_name = "hs071_SNOPT_setDVFromHist"
        optOptions = {
            "Print file": "{}.out".format(test_name),
            "Summary file": "{}_summary.out".format(test_name),
        }
        histFileName = "{}.hst".format(test_name)
        optProb = self.setup_optProb(xScale=1.5, conScale=1.2, objScale=32, offset=1.5)
        sol = self.optimize(optProb, "SNOPT", storeHistory=histFileName, optOptions=optOptions)
        self.check_solution(sol, self.tol["SNOPT"])
        hist = History(histFileName, flag="r")
        first = hist.getValues(names="xvars", callCounters="last", scale=False)
        x_final = first["xvars"][0]
        optProb = self.setup_optProb(xScale=0.5, conScale=4.8, objScale=0.1, offset=1.5)
        sol = self.optimize(optProb, "SNOPT", setDV=histFileName, storeHistory=histFileName, optOptions=optOptions)
        self.check_solution(sol, self.tol["SNOPT"])
        # Verify the history file
        hist = History(histFileName, flag="r")
        second = hist.getValues(names="xvars", scale=False)
        x_init = second["xvars"][0]
        assert_allclose(x_init, x_final, atol=1e-5, rtol=1e-5)
        # assert that this only took one major iteration
        # since we restarted from the optimum
        self.assertEqual(second["xvars"].shape, (1, 4))

    def test_slsqp_scaling_offset_optProb(self):
        """
        Test that scaling and offset works as expected
        Also test optProb stored in the history file is correct
        """
        test_name = "hs071_SLSQP_scaling_offset"
        histFileName = "{}.hst".format(test_name)
        optOptions = {"IFILE": "{}.out".format(test_name)}
        objScale = 4.2
        xScale = [2, 3, 4, 5]
        conScale = [0.6, 1.7]
        offset = [1, -2, 40, 2.5]
        optProb = self.setup_optProb(objScale=objScale, xScale=xScale, conScale=conScale, offset=offset)
        sol = self.optimize(optProb, "SLSQP", storeHistory=histFileName, optOptions=optOptions)
        self.check_solution(sol, self.tol["SLSQP"])
        # now we retrieve the history file, and check the scale=True option is indeed
        # scaling things correctly
        hist = History(histFileName, flag="r")
        orig_values = hist.getValues(callCounters="0", scale=False)
        optProb = hist.getOptProb()

        # check that the scales are stored properly
        for i, var in enumerate(optProb.variables["xvars"]):
            assert_allclose(xScale[i], var.scale, atol=1e-12, rtol=1e-12)
            assert_allclose(offset[i], var.offset, atol=1e-12, rtol=1e-12)
        for con in optProb.constraints:
            assert_allclose(conScale, optProb.constraints[con].scale, atol=1e-12, rtol=1e-12)
        for obj in optProb.objectives:
            assert_allclose(objScale, optProb.objectives[obj].scale, atol=1e-12, rtol=1e-12)

        # verify the scale option in getValues
        scaled_values = hist.getValues(callCounters="0", scale=True, stack=False)
        x = orig_values["xvars"][0]
        x_scaled = scaled_values["xvars"][0]
        assert_allclose(x_scaled, (x - offset) * xScale, atol=1e-12, rtol=1e-12)

        # now do the same but with stack=True
        stacked_values = hist.getValues(callCounters="0", scale=True, stack=True)
        x_scaled = stacked_values["xuser"][0]
        assert_allclose(x_scaled, (x - offset) * xScale, atol=1e-12, rtol=1e-12)

        # now we test objective and constraint scaling in getValues
        obj_orig = orig_values["obj"][0]
        obj_scaled = scaled_values["obj"][0]
        assert_allclose(obj_scaled, obj_orig * objScale, atol=1e-12, rtol=1e-12)
        con_orig = orig_values["con"][0]
        con_scaled = scaled_values["con"][0]
        assert_allclose(con_scaled, con_orig * conScale, atol=1e-12, rtol=1e-12)

    def test_snopt_informs(self):
        optProb = self.setup_optProb()
        test_name = "hs071_SNOPT"
        optOptions = {
            "Print file": "{}.out".format(test_name),
            "Summary file": "{}_summary.out".format(test_name),
        }
        optOptions["Major iterations limit"] = 1
        sol = self.optimize(optProb, "SNOPT", optOptions=optOptions)
        self.assertEqual(sol.optInform["value"], 32)

    def test_slsqp_informs(self):
        optProb = self.setup_optProb()
        optOptions = {"IFILE": "hs071_SLSQP.out"}
        # now we set max iteration to 1 and verify that we get a different inform
        optOptions["MAXIT"] = 1
        sol = self.optimize(optProb, "slsqp", optOptions=optOptions)
        self.assertEqual(sol.optInform["value"], 9)

    def test_nlpqlp_informs(self):
        optProb = self.setup_optProb()
        optOptions = {"iFile": "hs071_NLPQLP.out"}
        optOptions["maxIt"] = 1
        sol = self.optimize(optProb, "nlpqlp", optOptions=optOptions)
        self.assertEqual(sol.optInform["value"], 1)

    def test_ipopt_informs(self):
        optProb = self.setup_optProb()
        optOptions = {"print_level": 5, "output_file": "hs071_IPOPT.out"}
        # Test that the inform is -1 when iterations are too limited.
        optOptions["max_iter"] = 1
        sol = self.optimize(optProb, "ipopt", optOptions=optOptions)
        self.assertEqual(sol.optInform["value"], -1)
        # Test that the inform is -4 when max_cpu_time are too limited.
        optOptions["max_iter"] = 100
        optOptions["max_cpu_time"] = 0.001
        sol = self.optimize(optProb, "ipopt", optOptions=optOptions)
        self.assertEqual(sol.optInform["value"], -4)

    def test_psqp_informs(self):
        optProb = self.setup_optProb()
        optOptions = {"IFILE": "hs071_PSQP.out"}
        optOptions["MIT"] = 1
        sol = self.optimize(optProb, "psqp", optOptions=optOptions)
        self.assertEqual(sol.optInform["value"], 11)

    @parameterized.expand(
        [
            "SNOPT",
            "IPOPT",
            "SLSQP",
            "PSQP",
            "CONMIN",
            "NLPQLP",
            "ParOpt",
        ]
    )
    def test_optimization(self, optName):
        if optName in self.optOptions:
            optOptions = self.optOptions[optName]
        else:
            optOptions = {}
        optOptions = self.updateOptOptions(optName, optOptions)
        optProb = self.setup_optProb()
        sol = self.optimize(optProb, optName, optOptions=optOptions)
        # Check Solution
        self.check_solution(sol, self.tol[optName])
        # Check informs
        self.check_inform(sol, optName)


if __name__ == "__main__":
    unittest.main()
