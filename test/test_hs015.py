"""Test solution of problem HS15 from the Hock & Schittkowski collection"""

import unittest
import numpy as np
from numpy.testing import assert_allclose
from pyoptsparse import Optimization, OPT, History
from pyoptsparse.pyOpt_error import Error


class TestHS15(unittest.TestCase):

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

    def objfunc(self, xdict):
        self.nf += 1
        x = xdict["xvars"]
        funcs = {}
        funcs["obj"] = [100 * (x[1] - x[0] ** 2) ** 2 + (1 - x[0]) ** 2]
        conval = np.zeros(2, "D")
        conval[0] = x[0] * x[1]
        conval[1] = x[0] + x[1] ** 2
        funcs["con"] = conval
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

    def optimize(self, optName, tol, optOptions={}, storeHistory=False, hotStart=None):
        self.nf = 0  # number of function evaluations
        self.ng = 0  # number of gradient evaluations
        # Optimization Object
        optProb = Optimization("HS15 Constraint Problem", self.objfunc)

        # Design Variables
        lower = [-5.0, -5.0]
        upper = [0.5, 5.0]
        value = [-2, 1.0]
        optProb.addVarGroup("xvars", 2, lower=lower, upper=upper, value=value)

        # Constraints
        lower = [1.0, 0.0]
        upper = [None, None]
        optProb.addConGroup("con", 2, lower=lower, upper=upper)

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
                self.histFileName = "%s_hs015_Hist.hst" % (optName.lower())
            elif isinstance(storeHistory, str):
                self.histFileName = storeHistory
        else:
            self.histFileName = None

        sol = opt(optProb, sens=self.sens, storeHistory=self.histFileName, hotStart=hotStart)

        # Test printing solution to screen
        print(sol)

        # Check Solution
        self.fStar1 = 306.5
        self.fStar2 = 360.379767

        self.xStar1 = (0.5, 2.0)
        self.xStar2 = (-0.79212322, -1.26242985)

        dv = sol.getDVs()
        sol_xvars = [sol.variables["xvars"][i].value for i in range(2)]
        assert_allclose(sol_xvars, dv["xvars"], atol=tol, rtol=tol)
        # we check either optimum via try/except
        try:
            assert_allclose(sol.objectives["obj"].value, self.fStar1, atol=tol, rtol=tol)
            assert_allclose(dv["xvars"], self.xStar1, atol=tol, rtol=tol)
        except AssertionError:
            assert_allclose(sol.objectives["obj"].value, self.fStar2, atol=tol, rtol=tol)
            assert_allclose(dv["xvars"], self.xStar2, atol=tol, rtol=tol)

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
        self.assertEqual(hist.getConNames(), ["con"])
        self.assertEqual(hist.getObjNames(), ["obj"])
        dvInfo = hist.getDVInfo()
        self.assertEqual(len(dvInfo), 1)
        self.assertEqual(dvInfo["xvars"], hist.getDVInfo(key="xvars"))
        conInfo = hist.getConInfo()
        self.assertEqual(len(conInfo), 1)
        self.assertEqual(conInfo["con"], hist.getConInfo(key="con"))
        objInfo = hist.getObjInfo()
        self.assertEqual(len(objInfo), 1)
        self.assertEqual(objInfo["obj"], hist.getObjInfo(key="obj"))
        for key in ["lower", "upper", "scale"]:
            self.assertIn(key, dvInfo["xvars"])
            self.assertIn(key, conInfo["con"])
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
        test_name = "hs015_SNOPT"
        store_vars = ["step", "merit", "feasibility", "optimality", "penalty", "Hessian", "condZHZ", "slack", "lambda"]
        optOptions = {
            "Save major iteration variables": store_vars,
            "Print file": "{}.out".format(test_name),
            "Summary file": "{}_summary.out".format(test_name),
        }
        self.optimize_with_hotstart("SNOPT", 1e-12, optOptions=optOptions)

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

    def test_slsqp(self):
        optOptions = {"IFILE": "hs015_SLSQP.out"}
        self.optimize_with_hotstart("SLSQP", 1e-8, optOptions=optOptions)

    def test_nlpqlp(self):
        optOptions = {"iFile": "hs015_NLPQLP.out"}
        self.optimize_with_hotstart("NLPQLP", 1e-12, optOptions=optOptions)

    def test_ipopt(self):
        optOptions = {"output_file": "hs015_IPOPT.out"}
        self.optimize_with_hotstart("IPOPT", 1e-4, optOptions=optOptions)

    def test_paropt(self):
        optOptions = {"output_file": "hs015_ParOpt.out"}
        self.optimize_with_hotstart("ParOpt", 1e-6, optOptions=optOptions)

    def test_conmin(self):
        optOptions = {
            "DELFUN": 1e-10,
            "DABFUN": 1e-10,
            "IFILE": "hs015_CONMIN.out",
        }
        self.optimize_with_hotstart("CONMIN", 1e-10, optOptions=optOptions)

    def test_psqp(self):
        optOptions = {"IFILE": "hs015_PSQP.out"}
        self.optimize_with_hotstart("PSQP", 1e-12, optOptions=optOptions)


if __name__ == "__main__":
    unittest.main()
