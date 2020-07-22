"""Test solution of problem HS71 from the Hock & Schittkowski collection"""

import unittest
import numpy as np
from numpy.testing import assert_allclose
from pyoptsparse import Optimization, OPT, History
from pyoptsparse.pyOpt_error import Error


class TestHS71(unittest.TestCase):
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

    def optimize(
        self,
        optName,
        tol,
        optOptions={},
        storeHistory=False,
        setDV=None,
        xScale=1.0,
        objScale=1.0,
        conScale=1.0,
        offset=0.0,
        check_solution=True,
    ):
        # Optimization Object
        optProb = Optimization("HS071 Constraint Problem", self.objfunc)

        # Design Variables
        x0 = [1.0, 5.0, 5.0, 1.0]
        optProb.addVarGroup("xvars", 4, lower=1, upper=5, value=x0, scale=xScale, offset=offset)

        # Constraints
        optProb.addConGroup("con", 2, lower=[25, 40], upper=[None, 40], scale=conScale)

        # Objective
        optProb.addObj("obj", scale=objScale)

        # Optimizer
        try:
            opt = OPT(optName, options=optOptions)
        except Error:
            raise unittest.SkipTest("Optimizer not available:", optName)

        if isinstance(setDV, str):
            optProb.setDVsFromHistory(setDV)
        elif isinstance(setDV, dict):
            optProb.setDVs(setDV)
            outDV = optProb.getDVs()
            assert_allclose(setDV["xvars"], outDV["xvars"])

        sol = opt(optProb, sens=self.sens, storeHistory=storeHistory)

        # Check Solution
        if check_solution:
            self.fStar = 17.0140172
            self.xStar = (1.0, 4.743, 3.82115, 1.37941)
            self.lambdaStar = (0.55229366, -0.16146857)
            assert_allclose(sol.objectives["obj"].value, self.fStar, atol=tol, rtol=tol)
            assert_allclose(sol.xStar["xvars"], self.xStar, atol=tol, rtol=tol)

            if hasattr(sol, "lambdaStar"):
                assert_allclose(sol.lambdaStar["con"], self.lambdaStar, atol=tol, rtol=tol)
        return sol

    def test_slsqp_setDV(self):
        """
        Test that setDV works as expected, even with scaling/offset
        """
        test_name = "hs071_SLSQP_setDV"
        histFileName = "{}.hst".format(test_name)
        newDV = {"xvars": np.array([1, 4, 4, 1])}
        optOptions = {"IFILE": "{}.out".format(test_name)}
        self.optimize(
            "SLSQP",
            1e-5,
            xScale=1.5,
            conScale=1.2,
            objScale=32,
            offset=1.5,
            setDV=newDV,
            storeHistory=histFileName,
            optOptions=optOptions,
        )
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
        self.optimize(
            "snopt",
            1e-6,
            xScale=1.5,
            conScale=1.2,
            objScale=32,
            offset=1.5,
            storeHistory=histFileName,
            optOptions=optOptions,
        )
        hist = History(histFileName, flag="r")
        first = hist.getValues(names="xvars", callCounters="last", scale=False)
        x_final = first["xvars"][0]
        self.optimize(
            "snopt",
            1e-6,
            xScale=0.5,
            conScale=4.8,
            objScale=0.1,
            offset=1.5,
            setDV=histFileName,
            storeHistory=histFileName,
            optOptions=optOptions,
        )
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
        self.optimize(
            "slsqp",
            1e-6,
            objScale=objScale,
            xScale=xScale,
            conScale=conScale,
            storeHistory=histFileName,
            offset=offset,
            optOptions=optOptions,
        )

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

    def test_snopt(self):
        test_name = "hs071_SNOPT"
        optOptions = {
            "Print file": "{}.out".format(test_name),
            "Summary file": "{}_summary.out".format(test_name),
        }
        sol = self.optimize("snopt", 1e-6, optOptions=optOptions)
        self.assertEqual(sol.optInform["value"], 1)
        optOptions["Major iterations limit"] = 1
        sol = self.optimize("snopt", 1e-6, optOptions=optOptions, check_solution=False)
        self.assertEqual(sol.optInform["value"], 32)

    def test_slsqp(self):
        optOptions = {"IFILE": "hs071_SLSQP.out"}
        sol = self.optimize("slsqp", 1e-6, optOptions=optOptions)
        self.assertEqual(sol.optInform["value"], 0)
        # now we set max iteration to 1 and verify that we get a different inform
        optOptions["MAXIT"] = 1
        sol = self.optimize("slsqp", 1e-6, optOptions=optOptions, check_solution=False)
        self.assertEqual(sol.optInform["value"], 9)

    def test_nlpqlp(self):
        optOptions = {"iFile": "hs071_NLPQLP.out"}
        sol = self.optimize("nlpqlp", 1e-6, optOptions=optOptions)
        self.assertEqual(sol.optInform["value"], 0)
        optOptions["maxIt"] = 1
        sol = self.optimize("nlpqlp", 1e-6, optOptions=optOptions, check_solution=False)
        self.assertEqual(sol.optInform["value"], 1)

    def test_ipopt(self):
        optOptions = {"print_level": 5, "output_file": "hs071_IPOPT.out"}
        sol = self.optimize("ipopt", 1e-6, optOptions=optOptions)
        self.assertEqual(sol.optInform["value"], 0)
        self.assertEqual(sol.optInform["text"], "Solve Succeeded")
        # Test that the inform is -1 when iterations are too limited.
        optOptions["max_iter"] = 1
        sol = self.optimize("ipopt", 1e-6, optOptions=optOptions, check_solution=False)
        self.assertEqual(sol.optInform["value"], -1)
        self.assertEqual(sol.optInform["text"], "Maximum Iterations Exceeded")
        # Test that the inform is -4 when max_cpu_time are too limited.
        optOptions["max_iter"] = 100
        optOptions["max_cpu_time"] = 0.001
        sol = self.optimize("ipopt", 1e-6, optOptions=optOptions, check_solution=False)
        self.assertEqual(sol.optInform["value"], -4)
        self.assertEqual(sol.optInform["text"], "Maximum CpuTime Exceeded")

    def test_conmin(self):
        optOptions = {
            "DELFUN": 1e-9,
            "DABFUN": 1e-9,
            "IFILE": "hs071_CONMIN.out",
        }
        self.optimize("conmin", 1e-2, optOptions=optOptions)

    def test_psqp(self):
        optOptions = {"IFILE": "hs071_PSQP.out"}
        sol = self.optimize("psqp", 1e-6, optOptions=optOptions)
        self.assertEqual(sol.optInform["value"], 4)
        optOptions["MIT"] = 1
        sol = self.optimize("psqp", 1e-6, optOptions=optOptions, check_solution=False)
        self.assertEqual(sol.optInform["value"], 11)

    def test_paropt(self):
        optOptions = {"output_file": "hs071_ParOpt.out"}
        self.optimize("paropt", 1e-6, optOptions=optOptions)


if __name__ == "__main__":
    unittest.main()
