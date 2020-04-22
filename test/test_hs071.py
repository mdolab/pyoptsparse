"""Test solution of problem HS71 from the Hock & Schittkowski collection"""

import unittest
import numpy
from numpy.testing import assert_allclose
from pyoptsparse import Optimization, OPT, History

class TestHS71(unittest.TestCase):
    def objfunc(self, xdict):
        x = xdict['xvars']
        funcs = {}
        funcs['obj'] = x[0] * x[3] * (x[0] + x[1] + x[2]) + x[2]
        funcs['con'] = [x[0] * x[1] * x[2] * x[3],
                    x[0]*x[0] + x[1]*x[1] + x[2]*x[2] + x[3]*x[3]]
        fail = False
        return funcs, fail


    def sens(self, xdict, funcs):
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

    def optimize(self, optName, tol, optOptions={}, storeHistory=False, xScale=1.0, 
                 objScale=1.0, conScale=1.0, check_solution=True):
        # Optimization Object
        optProb = Optimization('HS071 Constraint Problem', self.objfunc)

        # Design Variables
        x0 = [1.0, 5.0, 5.0, 1.0]
        optProb.addVarGroup('xvars', 4, lower=1, upper=5, value=x0, scale=xScale)

        # Constraints
        optProb.addConGroup('con', 2, lower=[25, 40], upper=[None, 40], scale=conScale)

        # Objective
        optProb.addObj('obj', scale=objScale)

        # Optimizer
        try:
            opt = OPT(optName, options=optOptions)
        except:
            raise unittest.SkipTest('Optimizer not available:', optName)

        sol = opt(optProb, sens=self.sens, storeHistory=storeHistory)

        # Check Solution
        if check_solution:
            self.fStar = 17.0140172
            self.xStar = (1.0, 4.743, 3.82115, 1.37941)
            self.lambdaStar = (0.55229366, -0.16146857)
            assert_allclose(sol.objectives['obj'].value, self.fStar, atol = tol, rtol = tol)
            assert_allclose(sol.xStar['xvars'], self.xStar, atol = tol, rtol = tol)

            if hasattr(sol, 'lambdaStar'):
                assert_allclose(sol.lambdaStar['con'], self.lambdaStar, atol = tol, rtol = tol)
        return sol

    def test_snopt(self):
        self.optimize('snopt', 1E-6)
    
    def test_snopt_scaling(self):
        histFileName = 'snopt_scale_test.hst'
        objScale = 4.2
        xScale = [2,3,4,5]
        conScale = [0.6, 1.7]
        self.optimize('snopt', 1E-6, objScale=objScale, xScale=xScale, conScale=conScale, storeHistory=histFileName)

        # now we retrieve the history file, and check the scale=True option is indeed
        # scaling things correctly
        hist = History(histFileName, flag='r')
        last = hist.getValues(callCounters='last', scale=False)
        last_scaled = hist.getValues(callCounters='last', scale=True)
        x = last['xvars'][0]
        x_scaled = last_scaled['xvars'][0]
        assert_allclose(x_scaled/x,xScale, atol=1E-12, rtol=1E-12)

    def test_slsqp(self):
        self.optimize('slsqp', 1E-6)

    def test_nlpqlp(self):
        self.optimize('nlpqlp', 1E-6)

    def test_ipopt(self):
        opts = {}
        opts['print_level'] = 5
        sol = self.optimize('ipopt', 1E-6, optOptions=opts)
        self.assertEqual(sol.optInform['value'], 0)
        self.assertEqual(sol.optInform['text'], 'Solve Succeeded')
        # Test that the inform is -1 when iterations are too limited.
        opts['max_iter'] = 1
        sol = self.optimize('ipopt', 1E-6, optOptions=opts, check_solution=False)
        self.assertEqual(sol.optInform['value'], -1)
        self.assertEqual(sol.optInform['text'], 'Maximum Iterations Exceeded')
        # Test that the inform is -4 when max_cpu_time are too limited.
        opts['max_iter'] = 100
        opts['max_cpu_time'] = 0.001
        sol = self.optimize('ipopt', 1E-6, optOptions=opts, check_solution=False)
        self.assertEqual(sol.optInform['value'], -4)
        self.assertEqual(sol.optInform['text'], 'Maximum CpuTime Exceeded')

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
