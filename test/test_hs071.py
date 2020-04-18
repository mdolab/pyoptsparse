"""Test solution of problem HS71 from the Hock & Schittkowski collection"""
from __future__ import print_function

import unittest

import numpy
from pyoptsparse import Optimization, OPT

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

    def optimize(self, optName, optOptions={}, storeHistory=False, places=5, xScale=1.0, objScale=1.0, conScale=1.0):
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

        sol = opt(optProb, sens=self.sens)

        # Check Solution
        self.fStar = 17.0140172
        self.xStar = (1.0, 4.743, 3.82115, 1.37941)
        self.lambdaStar = (0.55229366, -0.16146857)
        self.assertAlmostEqual(sol.objectives['obj'].value, self.fStar, places=places)
        for i in range(4):
            self.assertAlmostEqual(sol.xStar['xvars'][i], self.xStar[i], places=places)

        if hasattr(sol, 'lambdaStar'):
            for i in range(2):
                self.assertAlmostEqual(sol.lambdaStar['con'][i], self.lambdaStar[i], places=places)

    def test_snopt(self):
        self.optimize('snopt')
    
    def test_snopt_scaling(self):
        self.optimize('snopt', objScale=4.2, xScale=[2,3,4,5], conScale=[0.6, 1.7])
        

    def test_slsqp(self):
        self.optimize('slsqp')

    def test_nlpqlp(self):
        self.optimize('nlpqlp')

    def test_ipopt(self):
        self.optimize('ipopt')

    def test_conmin(self):
        opts = {'DELFUN' : 1e-9,
                'DABFUN' : 1e-9}
        self.optimize('conmin', optOptions=opts, places=2)

    def test_psqp(self):
        self.optimize('psqp')

    def test_paropt(self):
        self.optimize('paropt')

if __name__ == "__main__":
    unittest.main()
