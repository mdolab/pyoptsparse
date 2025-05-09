"""Test solution of problem HS15 from the Hock & Schittkowski collection"""

# Standard Python modules
import os
import unittest

# External modules
from baseclasses.utils import readPickle, writePickle
import numpy as np
from parameterized import parameterized

# First party modules
from pyoptsparse import OPT, History, Optimization

# Local modules
from testing_utils import OptTest


class TestHS15(OptTest):
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

    name = "HS015"
    DVs = {"xvars"}
    cons = {"con"}
    objs = {"obj"}
    extras = {"extra1", "extra2"}
    fStar = [
        306.5,
        360.379767,
    ]
    xStar = [
        {"xvars": (0.5, 2.0)},
        {"xvars": (-0.79212322, -1.26242985)},
    ]
    tol = {
        "SLSQP": 1e-5,
        "NLPQLP": 1e-12,
        "IPOPT": 1e-4,
        "ParOpt": 1e-6,
        "CONMIN": 1e-10,
        "PSQP": 5e-12,
    }
    optOptions = {}

    def objfunc(self, xdict):
        self.nf += 1
        x = xdict["xvars"]
        funcs = {}
        funcs["obj"] = [100 * (x[1] - x[0] ** 2) ** 2 + (1 - x[0]) ** 2]
        conval = np.zeros(2, "D")
        conval[0] = x[0] * x[1]
        conval[1] = x[0] + x[1] ** 2
        funcs["con"] = conval
        # extra keys
        funcs["extra1"] = 0.0
        funcs["extra2"] = 1.0
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

    def setup_optProb(self):
        # Optimization Object
        self.optProb = Optimization("HS15 Constraint Problem", self.objfunc)

        # Design Variables
        lower = [-5.0, -5.0]
        upper = [0.5, 5.0]
        value = [-2, 1.0]
        self.optProb.addVarGroup("xvars", 2, lower=lower, upper=upper, value=value)

        # Constraints
        lower = [1.0, 0.0]
        upper = [None, None]
        self.optProb.addConGroup("con", 2, lower=lower, upper=upper)

        # Objective
        self.optProb.addObj("obj")

    def test_snopt(self):
        self.optName = "SNOPT"
        self.setup_optProb()
        store_vars = ["Hessian", "slack", "lambda", "penalty_vector", "nS", "BSwap", "maxVi"]
        optOptions = {"Save major iteration variables": store_vars}
        self.optimize_with_hotstart(1e-12, optOptions=optOptions)

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
        # dv = sol.getDVs()
        # sol_xvars = [sol.variables["xvars"][i].value for i in range(2)]
        # assert_allclose(sol_xvars, dv["xvars"], atol=tol, rtol=tol)

    @parameterized.expand(["SLSQP", "PSQP", "CONMIN", "NLPQLP", "ParOpt"])
    def test_optimization(self, optName):
        self.optName = optName
        self.setup_optProb()
        optOptions = self.optOptions.pop(optName, None)
        sol = self.optimize(optOptions=optOptions)
        # Check Solution
        self.assert_solution_allclose(sol, self.tol[optName])
        # Check informs
        self.assert_inform_equal(sol)

    def test_ipopt(self):
        self.optName = "IPOPT"
        self.setup_optProb()
        optOptions = self.optOptions.pop(self.optName, None)
        sol = self.optimize(optOptions=optOptions, storeHistory=True)
        # Check Solution
        self.assert_solution_allclose(sol, self.tol[self.optName])
        # Check informs
        self.assert_inform_equal(sol)

        # Check iteration counters
        hist = History(self.histFileName, flag="r")
        data_init = hist.read(0)
        self.assertEqual(0, data_init["iter"])
        data_last = hist.read(hist.read("last"))
        self.assertGreater(data_last["iter"], 0)

        # Make sure there is no duplication in objective history
        data = hist.getValues(names=["obj"])
        objhis_len = data["obj"].shape[0]
        for i in range(objhis_len - 1):
            self.assertNotEqual(data["obj"][i], data["obj"][i + 1])

    def test_snopt_hotstart(self):
        self.optName = "SNOPT"
        self.setup_optProb()
        sol, restartDict = self.optimize(optOptions={"Return work arrays": True})
        # Check Solution
        self.assert_solution_allclose(sol, 1e-12)
        # Check informs
        self.assert_inform_equal(sol)
        # Check restartDict
        self.assertEqual({"cw", "iw", "rw", "xs", "hs", "pi"}, set(restartDict.keys()))

        # Now optimize again, but using the hotstart
        self.setup_optProb()
        self.nf = 0
        self.ng = 0
        opt = OPT(self.optName, options={"Start": "Hot", "Verify level": -1})
        sol = opt(self.optProb, sens=self.sens, restartDict=restartDict)
        # Check Solution
        self.assert_solution_allclose(sol, 1e-12)
        # Should only take one major iteration
        self.assertEqual(self.nf, 1)
        self.assertEqual(self.ng, 1)

    @staticmethod
    def my_snstop(iterDict):
        """manually terminate SNOPT after 1 major iteration"""
        if iterDict["nMajor"] == 1:
            return 1
        return 0

    def test_snopt_snstop(self):
        self.optName = "SNOPT"
        self.setup_optProb()
        optOptions = {
            "snSTOP function handle": self.my_snstop,
        }
        sol = self.optimize(optOptions=optOptions, storeHistory=True)
        # Check informs
        # we should get 70/74
        self.assert_inform_equal(sol, optInform=74)

    def test_snopt_snstop_restart(self):
        pickleFile = "restart.pickle"

        def my_snstop_restart(iterDict, restartDict):
            # Save the restart dictionary
            writePickle(pickleFile, restartDict)

            # Exit after 5 major iterations
            if iterDict["nMajor"] == 5:
                return 1

            return 0

        # Run the optimization for 5 major iterations
        self.optName = "SNOPT"
        self.setup_optProb()
        optOptions = {
            "snSTOP function handle": my_snstop_restart,
            "snSTOP arguments": ["restartDict"],
        }
        sol = self.optimize(optOptions=optOptions, storeHistory=True)

        # Check that the optimization exited with 74
        self.assert_inform_equal(sol, optInform=74)

        # Read the restart dictionary pickle file saved by snstop
        restartDict = readPickle(pickleFile)

        # Now optimize again but using the restart dictionary
        self.setup_optProb()
        opt = OPT(
            self.optName,
            options={
                "Start": "Hot",
                "Verify level": -1,
                "snSTOP function handle": my_snstop_restart,
                "snSTOP arguments": ["restartDict"],
            },
        )
        histFile = "restart.hst"
        sol = opt(self.optProb, sens=self.sens, storeHistory=histFile, restartDict=restartDict)

        # Check that the optimization converged in fewer than 5 more major iterations
        self.assert_solution_allclose(sol, 1e-12)
        self.assert_inform_equal(sol, optInform=1)

        # Delete the pickle and history files
        os.remove(pickleFile)
        os.remove(histFile)

    def test_snopt_work_arrays_save(self):
        # Run the optimization for 5 major iterations
        self.optName = "SNOPT"
        self.setup_optProb()
        pickleFile = "work_arrays_save.pickle"
        optOptions = {
            "snSTOP function handle": self.my_snstop,
            "Work arrays save file": pickleFile,
        }
        sol = self.optimize(optOptions=optOptions, storeHistory=True)

        # Read the restart dictionary pickle file saved by snstop
        restartDict = readPickle(pickleFile)

        # Now optimize again but using the restart dictionary
        self.setup_optProb()
        opt = OPT(
            self.optName,
            options={
                "Start": "Hot",
                "Verify level": -1,
            },
        )
        histFile = "work_arrays_save.hst"
        sol = opt(self.optProb, sens=self.sens, storeHistory=histFile, restartDict=restartDict)

        # Check that the optimization converged
        self.assert_solution_allclose(sol, 1e-12)
        self.assert_inform_equal(sol, optInform=1)

        # Delete the pickle and history files
        os.remove(pickleFile)
        os.remove(histFile)

    @parameterized.expand(["SNOPT", "IPOPT"])
    def test_failed_initial_func(self, optName):
        # Test if we get the correct inform when the initial function call fails
        def failed_fun(x_dict):
            funcs = {"obj": 0.0, "con": [np.nan, np.nan]}
            fail = True
            return funcs, fail

        self.optName = optName
        self.setup_optProb()
        # swap obj to report NaN
        self.optProb.objFun = failed_fun
        sol = self.optimize(optOptions={}, storeHistory=True)
        if self.optName == "SNOPT":
            inform_ref = 61
        elif self.optName == "IPOPT":
            inform_ref = -13
        self.assert_inform_equal(sol, optInform=inform_ref)
        # make sure empty history does not error out
        hist = History(self.histFileName, flag="r")
        hist.getValues()

    @parameterized.expand(["SNOPT", "IPOPT"])
    def test_failed_initial_sens(self, optName):
        # Test if we get the correct inform when the initial sensitivity call fails
        def failed_sens(xdict, funcs):
            funcsSens = {}
            funcsSens["obj"] = {"xvars": [np.nan, np.nan]}
            funcsSens["con"] = {"xvars": [[np.nan, np.nan], [np.nan, np.nan]]}
            fail = True
            return funcsSens, fail

        self.optName = optName
        self.setup_optProb()
        sol = self.optimize(sens=failed_sens, optOptions={}, storeHistory=True)
        if self.optName == "SNOPT":
            inform_ref = 61
        elif self.optName == "IPOPT":
            inform_ref = -13
        self.assert_inform_equal(sol, optInform=inform_ref)
        # make sure empty history does not error out
        hist = History(self.histFileName, flag="r")
        hist.getValues()

    @parameterized.expand(["SNOPT", "IPOPT"])
    def test_func_eval_failure(self, optName):
        # Test if optimizer back-tracks after function evaluation failure and keeps going

        # This function is the same as original except it will fail at 3rd call
        def objFun_eval_failure(x_dict):
            self.nf += 1
            x = x_dict["xvars"]
            funcs = {}
            funcs["obj"] = [100 * (x[1] - x[0] ** 2) ** 2 + (1 - x[0]) ** 2]
            conval = np.zeros(2, "D")
            conval[0] = x[0] * x[1]
            conval[1] = x[0] + x[1] ** 2
            funcs["con"] = conval
            # extra keys
            funcs["extra1"] = 0.0
            funcs["extra2"] = 1.0
            fail = False
            if self.nf == 3:
                funcs["obj"] = [np.nan]
                funcs["con"] = [np.nan, np.nan]
                fail = True
            return funcs, fail

        self.optName = optName
        self.setup_optProb()
        # swap obj to simulate eval failure
        self.optProb.objFun = objFun_eval_failure
        sol = self.optimize(optOptions={}, storeHistory=True)
        # check solution and informs
        self.assert_solution_allclose(sol, 1e-5)
        self.assert_inform_equal(sol)


if __name__ == "__main__":
    unittest.main()
