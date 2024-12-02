"""Test solution of problem HS71 from the Hock & Schittkowski collection"""

# Standard Python modules
import unittest

# External modules
import numpy as np
from numpy.testing import assert_allclose
from parameterized import parameterized

# First party modules
from pyoptsparse import History, Optimization

# Local modules
from testing_utils import OptTest


class TestHS71(OptTest):
    # Optimization problem definition
    name = "hs071"
    DVs = {"xvars"}
    cons = {"con"}
    objs = {"obj"}
    fStar = 17.0140172
    xStar = {"xvars": (1.0, 4.743, 3.82115, 1.37941)}
    lambdaStar = {"con": (0.55229366, -0.16146857)}

    # Tolerances
    tol = {
        "SNOPT": 1e-6,
        "IPOPT": 1e-6,
        "NLPQLP": 1e-6,
        "SLSQP": 1e-6,
        "CONMIN": 1e-3,
        "PSQP": 1e-6,
        "ParOpt": 1e-6,
    }
    optOptions = {
        "CONMIN": {
            "DELFUN": 1e-10,
            "DABFUN": 1e-10,
        }
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

    def setup_optProb(self, xScale=1.0, objScale=1.0, conScale=1.0, offset=0.0):
        # Optimization Object
        self.optProb = Optimization("HS071 Constraint Problem", self.objfunc, sens=self.sens)

        # Design Variables
        x0 = [1.0, 5.0, 5.0, 1.0]
        self.optProb.addVarGroup("xvars", 4, lower=1, upper=5, value=x0, scale=xScale, offset=offset)

        # Constraints
        self.optProb.addConGroup("con", 2, lower=[25, 40], upper=[None, 40], scale=conScale)

        # Objective
        self.optProb.addObj("obj", scale=objScale)

    def test_slsqp_setDV(self):
        """
        Test that setDV works as expected, even with scaling/offset
        """
        self.optName = "SLSQP"
        histFileName = "hs071_SLSQP_setDV.hst"
        newDV = {"xvars": np.array([1, 4, 4, 1])}
        self.setup_optProb(xScale=1.5, conScale=1.2, objScale=32, offset=1.5)
        sol = self.optimize(setDV=newDV, storeHistory=histFileName)
        self.assert_solution_allclose(sol, self.tol["SLSQP"])
        # Verify the history file
        hist = History(histFileName, flag="r")
        init = hist.getValues(names="xvars", callCounters="0", scale=False)
        x_init = init["xvars"][0]
        assert_allclose(x_init, newDV["xvars"], atol=1e-5, rtol=1e-5)

    def test_snopt_setDVFromHist(self):
        """
        Test that setDVFromHistory works as expected, even with scaling/offset
        """
        self.optName = "SNOPT"
        histFileName = "hs071_SNOPT_setDVFromHist.hst"
        self.setup_optProb(xScale=1.5, conScale=1.2, objScale=32, offset=1.5)
        sol = self.optimize(storeHistory=histFileName)
        self.assert_solution_allclose(sol, self.tol["SNOPT"])
        hist = History(histFileName, flag="r")
        first = hist.getValues(names="xvars", callCounters="last", scale=False)
        x_final = first["xvars"][0]
        self.setup_optProb(xScale=0.5, conScale=4.8, objScale=0.1, offset=1.5)
        sol = self.optimize(setDV=histFileName, storeHistory=histFileName)
        self.assert_solution_allclose(sol, self.tol["SNOPT"])
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
        self.optName = "SLSQP"
        histFileName = "hs071_SLSQP_scaling_offset.hst"
        objScale = 4.2
        xScale = [2, 3, 4, 5]
        conScale = [0.6, 1.7]
        offset = [1, -2, 40, 2.5]
        self.setup_optProb(objScale=objScale, xScale=xScale, conScale=conScale, offset=offset)
        sol = self.optimize(storeHistory=histFileName)
        self.assert_solution_allclose(sol, self.tol["SLSQP"])
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
        self.optName = "SNOPT"
        self.setup_optProb()
        sol = self.optimize(optOptions={"Major iterations limit": 1})
        self.assert_inform_equal(sol, 32)

    def test_slsqp_informs(self):
        self.optName = "SLSQP"
        self.setup_optProb()
        # now we set max iteration to 1 and verify that we get a different inform
        sol = self.optimize(optOptions={"MAXIT": 1})
        self.assert_inform_equal(sol, 9)

    def test_nlpqlp_informs(self):
        self.optName = "NLPQLP"
        self.setup_optProb()
        sol = self.optimize(optOptions={"maxIt": 1})
        self.assert_inform_equal(sol, 1)

    def test_ipopt_informs(self):
        self.optName = "IPOPT"
        self.setup_optProb()
        # Test that the inform is -1 when iterations are too limited.
        sol = self.optimize(optOptions={"max_iter": 1})
        self.assert_inform_equal(sol, -1)
        # Test that the inform is -4 when max_cpu_time are too limited.
        sol = self.optimize(optOptions={"max_cpu_time": 0.001})
        self.assert_inform_equal(sol, -4)

    def test_psqp_informs(self):
        self.optName = "PSQP"
        self.setup_optProb()
        sol = self.optimize(optOptions={"MIT": 1})
        self.assert_inform_equal(sol, 11)

    @parameterized.expand(["SNOPT", "IPOPT", "SLSQP", "PSQP", "CONMIN", "NLPQLP", "ParOpt"])
    def test_optimization(self, optName):
        self.optName = optName
        self.setup_optProb()
        optOptions = self.optOptions.pop(optName, None)
        sol = self.optimize(optOptions=optOptions)
        # Check Solution
        lambda_sign = -1.0 if optName == "IPOPT" else 1.0
        self.assert_solution_allclose(sol, self.tol[optName], lambda_sign=lambda_sign)
        # Check informs
        self.assert_inform_equal(sol)
        # Check the lagrange multipliers in the solution text
        lines = str(sol).split("\n")
        constraint_header_line_num = [i for i, line in enumerate(lines) if "Constraints" in line][0]
        con1_line_num = constraint_header_line_num + 2
        con2_line_num = constraint_header_line_num + 3
        lambda_con1 = float(lines[con1_line_num].split()[-1])
        lambda_con2 = float(lines[con2_line_num].split()[-1])
        if optName in ("IPOPT", "SNOPT", "ParOpt"):
            # IPOPT returns Lagrange multipliers with opposite sign than SNOPT and ParOpt
            lambda_sign = -1.0 if optName == "IPOPT" else 1.0
            assert_allclose(
                [lambda_con1, lambda_con2],
                lambda_sign * np.asarray(self.lambdaStar[0]["con"]),
                rtol=1.0e-5,
                atol=1.0e-5,
            )
        else:
            assert_allclose([lambda_con1, lambda_con2], [9.0e100, 9.0e100], rtol=1.0e-5, atol=1.0e-5)


if __name__ == "__main__":
    unittest.main()
