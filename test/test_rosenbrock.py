"""Test solution of Rosenbrock problem"""

import unittest
import numpy as np
from numpy.testing import assert_allclose
from pyoptsparse import Optimization, OPT, History
from pyoptsparse.pyOpt_error import Error


class TestRosenbrock(unittest.TestCase):

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

    def optimize(self, optName, tol, optOptions={}, storeHistory=False, hotStart=None):
        self.nf = 0  # number of function evaluations
        self.ng = 0  # number of gradient evaluations
        # Optimization Object

        optProb = Optimization("Rosenbrock Problem", self.objfunc)

        n = 4  # Number of design variables
        np.random.seed(10)
        value = np.random.normal(size=n)

        lower = np.ones(n) * -50
        upper = np.ones(n) * 50
        optProb.addVarGroup("xvars", n, lower=lower, upper=upper, value=value)

        # Objective
        optProb.addObj("obj")

        # Check optimization problem:
        print(optProb)

        # Optimizer
        try:
            opt = OPT(optName, options=optOptions)
        except Error:
            raise unittest.SkipTest("Optimizer not available:", optName)

        # Solution
        if storeHistory is not None:
            if storeHistory is True:
                self.histFileName = "%s_Rsbrk_Hist.hst" % (optName.lower())
            elif isinstance(storeHistory, str):
                self.histFileName = storeHistory
        else:
            self.histFileName = None

        sol = opt(optProb, sens=self.sens, storeHistory=self.histFileName, hotStart=hotStart)

        # Test printing solution to screen
        print(sol)

        # Check Solution
        self.fStar1 = 0.0

        self.xStar1 = np.ones(n)

        dv = sol.getDVs()
        sol_xvars = [sol.variables["xvars"][i].value for i in range(n)]

        assert_allclose(sol_xvars, dv["xvars"], atol=tol, rtol=tol)
        # we check either optimum via try/except
        # try:
        assert_allclose(sol.objectives["obj"].value, self.fStar1, atol=tol, rtol=tol)
        assert_allclose(dv["xvars"], self.xStar1, atol=tol, rtol=tol)
        # except AssertionError:
        #     pass

    def check_hist_file(self, optimizer, tol):
        """
        We check the history file here along with the API
        """
        hist = History(self.histFileName, flag="r")
        # Metadata checks
        metadata = hist.getMetadata()
        self.assertEqual(metadata["optimizer"], optimizer)
        metadata_def_keys = ["optName", "optOptions", "nprocs", "startTime", "endTime", "optTime", "version"]
        for key in metadata_def_keys:
            self.assertIn(key, metadata)
        hist.getOptProb()

        # Info checks
        self.assertEqual(hist.getDVNames(), ["xvars"])
        # self.assertEqual(hist.getConNames(), ["con"])
        self.assertEqual(hist.getObjNames(), ["obj"])
        dvInfo = hist.getDVInfo()
        self.assertEqual(len(dvInfo), 1)
        self.assertEqual(dvInfo["xvars"], hist.getDVInfo(key="xvars"))
        conInfo = hist.getConInfo()
        self.assertEqual(len(conInfo), 0)
        # self.assertEqual(conInfo["con"], hist.getConInfo(key="con"))
        objInfo = hist.getObjInfo()
        self.assertEqual(len(objInfo), 1)
        self.assertEqual(objInfo["obj"], hist.getObjInfo(key="obj"))
        for key in ["lower", "upper", "scale"]:
            self.assertIn(key, dvInfo["xvars"])
            # self.assertIn(key, conInfo["con"])
        self.assertIn("scale", objInfo["obj"])

        # callCounter checks
        callCounters = hist.getCallCounters()
        last = hist.read("last")  # 'last' key should be present
        self.assertIn(last, callCounters)

        # iterKey checks
        iterKeys = hist.getIterKeys()
        for key in ["xuser", "fail", "isMajor"]:
            self.assertIn(key, iterKeys)

        # this check is only used for optimizers that guarantee '0' and 'last' contain funcs
        if optimizer in ["SNOPT", "SLSQP", "PSQP"]:
            val = hist.getValues(callCounters=["0", "last"], stack=True)
            self.assertEqual(val["isMajor"].size, 2)
            self.assertTrue(val["isMajor"][0])  # the first callCounter must be a major iteration
            self.assertTrue(val["isMajor"][-1])  # the last callCounter must be a major iteration
            # check optimum stored in history file against xstar
            assert_allclose(val["xuser"][-1], self.xStar1, atol=tol, rtol=tol)

    def optimize_with_hotstart(self, optName, tol, optOptions={}):
        """
        This code will perform 4 optimizations, one real opt and three restarts.
        In this process, it will check various combinations of storeHistory and hotStart filenames.
        It will also call `check_hist_file` after the first optimization.
        """

        self.optimize(optName, tol, storeHistory=True, optOptions=optOptions)
        self.assertGreater(self.nf, 0)
        self.assertGreater(self.ng, 0)
        self.check_hist_file(optName, tol)

        # re-optimize with hotstart
        self.optimize(optName, tol, storeHistory=False, hotStart=self.histFileName, optOptions=optOptions)
        # we should have zero actual function/gradient evaluations
        self.assertEqual(self.nf, 0)
        self.assertEqual(self.ng, 0)
        # another test with hotstart, this time with storeHistory = hotStart
        self.optimize(optName, tol, storeHistory=True, hotStart=self.histFileName, optOptions=optOptions)
        # we should have zero actual function/gradient evaluations
        self.assertEqual(self.nf, 0)
        self.assertEqual(self.ng, 0)
        # final test with hotstart, this time with a different storeHistory
        self.optimize(
            optName,
            tol,
            storeHistory="{}_new_hotstart.hst".format(optName),
            hotStart=self.histFileName,
            optOptions=optOptions,
        )
        # we should have zero actual function/gradient evaluations
        self.assertEqual(self.nf, 0)
        self.assertEqual(self.ng, 0)

    def test_snopt(self):
        test_name = "Rsbrk_SNOPT"
        store_vars = ["step", "merit", "feasibility", "optimality", "penalty", "Hessian", "condZHZ", "slack", "lambda"]
        optOptions = {
            "Save major iteration variables": store_vars,
            "Print file": "{}.out".format(test_name),
            "Summary file": "{}_summary.out".format(test_name),
        }
        self.optimize_with_hotstart("SNOPT", 1e-8, optOptions=optOptions)

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

    def test_slsqp(self):
        optOptions = {"ACC": 1e-10, "IFILE": "Rsbrk_SLSQP.out"}
        self.optimize_with_hotstart("SLSQP", 1e-6, optOptions=optOptions)

    # def test_nlpqlp(self):  ######### requires the same fix for gcon!!
    #     optOptions = {"iFile": "Rsbrk_NLPQLP.out"}
    #     self.optimize_with_hotstart("NLPQLP", 1e-8, optOptions=optOptions)

    def test_ipopt(self):
        optOptions = {"output_file": "Rsbrk_IPOPT.out"}
        self.optimize_with_hotstart("IPOPT", 1e-6, optOptions=optOptions)

    def test_paropt(self):
        optOptions = {"output_file": "Rsbrk_ParOpt.out"}
        self.optimize_with_hotstart("ParOpt", 1e-8, optOptions=optOptions)

    def test_conmin(self):
        optOptions = {
            "DELFUN": 1e-10,
            "DABFUN": 1e-10,
            "IFILE": "Rsbrk_CONMIN.out",
        }
        self.optimize_with_hotstart("CONMIN", 1e-9, optOptions=optOptions)

    def test_psqp(self):
        optOptions = {"IFILE": "Rsbrk_PSQP.out"}
        self.optimize_with_hotstart("PSQP", 1e-8, optOptions=optOptions)


if __name__ == "__main__":
    unittest.main()
