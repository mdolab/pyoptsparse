"""Test solution of problem HS71 from the Hock & Schittkowski collection"""
from __future__ import print_function

import unittest
import numpy
from numpy.testing import assert_allclose
from pyoptsparse import Optimization, OPT, History

class TestHS71(unittest.TestCase):
    def objfunc(self, xdict):
        self.nf += 1
        x = xdict['xvars']
        funcs = {}
        funcs['obj'] = x[0] * x[3] * (x[0] + x[1] + x[2]) + x[2]
        funcs['con'] = [x[0] * x[1] * x[2] * x[3],
                    x[0]*x[0] + x[1]*x[1] + x[2]*x[2] + x[3]*x[3]]
        fail = False
        return funcs, fail


    def sens(self, xdict, funcs):
        self.ng += 1
        x = xdict['xvars']
        funcsSens = {}
        funcsSens['obj'] = {'xvars': numpy.array(
                [x[0] * x[3] + x[3] * (x[0] + x[1] + x[2]) ,
                x[0] * x[3],
                x[0] * x[3] + 1.0,
                x[0] * (x[0] + x[1] + x[2])
                ])}
        jac = [[x[1]*x[2]*x[3], x[0]*x[2]*x[3], x[0]*x[1]*x[3], x[0]*x[1]*x[2]],
            [2.0*x[0], 2.0*x[1], 2.0*x[2], 2.0*x[3]]]
        funcsSens['con'] = {'xvars': jac}
        fail = False
        return funcsSens, fail

    def optimize(self, optName, tol, optOptions={}, storeHistory=False, setDV=None,
        xScale=1.0, objScale=1.0, conScale=1.0, offset=0.0):
        self.nf = 0 # number of function evaluations
        self.ng = 0 # number of gradient evaluations

        # Optimization Object
        optProb = Optimization('HS071 Constraint Problem', self.objfunc)

        # Design Variables
        x0 = [1.0, 5.0, 5.0, 1.0]
        optProb.addVarGroup('xvars', 4, lower=1, upper=5, value=x0, scale=xScale, offset=offset)

        # Constraints
        optProb.addConGroup('con', 2, lower=[25, 40], upper=[None, 40], scale=conScale)

        # Objective
        optProb.addObj('obj', scale=objScale)
        optProb.finalizeDesignVariables()

        # Optimizer
        try:
            opt = OPT(optName, options=optOptions)
        except:
            raise unittest.SkipTest('Optimizer not available:', optName)

        if isinstance(setDV, str):
            optProb.setDVsFromHistory(setDV)
        elif isinstance(setDV, dict):
            optProb.setDVs(setDV)

        sol = opt(optProb, sens=self.sens, storeHistory=storeHistory)

        # Check Solution
        self.fStar = 17.0140172
        self.xStar = (1.0, 4.743, 3.82115, 1.37941)
        self.lambdaStar = (0.55229366, -0.16146857)
        assert_allclose(sol.objectives['obj'].value, self.fStar, atol = tol, rtol = tol)
        assert_allclose(sol.xStar['xvars'], self.xStar, atol = tol, rtol = tol)

        if hasattr(sol, 'lambdaStar'):
            assert_allclose(sol.lambdaStar['con'], self.lambdaStar, atol = tol, rtol = tol)

    def test_snopt(self):
        self.optimize('snopt', 1E-6)

    def test_slsqp_setDV(self):
        histFileName = 'slsqp_test_DV.hst'
        self.optimize('slsqp', 1E-6, xScale=1.5, conScale=1.2, objScale=32, offset=1.5, storeHistory=histFileName)
        self.optimize('slsqp', 1E-6, xScale=0.5, conScale=4.8, objScale=0.1, offset=1.5, setDV=histFileName)
        # optimization should exit after one iteration
        self.assertEqual(self.ng, 1)
        self.assertEqual(self.nf, 1)

    def test_slsqp_scaling_offset(self):
        histFileName = 'slsqp_scale_offset.hst'
        objScale = 4.2
        xScale = [2,3,4,5]
        conScale = [0.6, 1.7]
        offset = [1,-2,40,2.5]
        self.optimize('slsqp', 1E-6, objScale=objScale, xScale=xScale, conScale=conScale, storeHistory=histFileName, offset=offset)

        # now we retrieve the history file, and check the scale=True option is indeed
        # scaling things correctly
        hist = History(histFileName, flag='r')
        last = hist.getValues(callCounters='0', scale=False)
        optProb = hist.getOptProb()

        # check that the scales are stored properly
        for i, var in enumerate(optProb.variables['xvars']):
            assert_allclose(xScale[i], var.scale, 1E-12)
            assert_allclose(offset[i], var.offset, 1E-12)
        for con in optProb.constraints:
            assert_allclose(conScale, optProb.constraints[con].scale, 1E-12)
        for obj in optProb.objectives:
            assert_allclose(objScale, optProb.objectives[obj].scale, 1E-12)

        # verify the scale option in getValues
        last_scaled = hist.getValues(callCounters='0', scale=True)
        x = last['xvars'][0]
        x_scaled = last_scaled['xvars'][0]
        assert_allclose(x_scaled, (x - offset)*xScale, atol=1E-12, rtol=1E-12)

    def test_slsqp(self):
        self.optimize('slsqp', 1E-6)

    def test_nlpqlp(self):
        self.optimize('nlpqlp', 1E-6)

    def test_ipopt(self):
        self.optimize('ipopt', 1E-6)

    def test_conmin(self):
        opts = {'DELFUN' : 1e-9,
                'DABFUN' : 1e-9}
        self.optimize('conmin', 1E-2, optOptions=opts)

    def test_psqp(self):
        self.optimize('psqp', 1E-6)

    def test_paropt(self):
        self.optimize('paropt', 1E-6)

if __name__ == "__main__":
    unittest.main()
