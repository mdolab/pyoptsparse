"""Test solution of problem HS71 from the Hock & Schittkowski collection"""
from __future__ import print_function

import unittest

import numpy
from pyoptsparse import Optimization, OPT


def objfunc(xdict):
    x = xdict['xvars']
    funcs = {}
    funcs['obj'] = x[0] * x[3] * (x[0] + x[1] + x[2]) + x[2]
    funcs['con'] = [x[0] * x[1] * x[2] * x[3],
                   x[0]*x[0] + x[1]*x[1] + x[2]*x[2] + x[3]*x[3]]
    fail = False
    return funcs, fail


def sens(xdict, funcs):
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


class TestHS71(unittest.TestCase):

    def optimize(self, optName, optOptions={}, storeHistory=False, places=5):
        # Optimization Object
        optProb = Optimization('HS071 Constraint Problem', objfunc)

        # Design Variables
        x0 = [1.0, 5.0, 5.0, 1.0]
        optProb.addVarGroup('xvars', 4, lower=1, upper=5, value=x0)

        # Constraints
        optProb.addConGroup('con', 2, lower=[25, 40], upper=[1e19, 40])

        # Objective
        optProb.addObj('obj')

        # Optimizer
        try:
            opt = OPT(optName, options=optOptions)
        except:
            raise unittest.SkipTest('Optimizer not available:', optName)

        sol = opt(optProb, sens=sens)

        # Check Solution
        self.assertAlmostEqual(sol.objectives['obj'].value, 17.0140172, places=places)

        self.assertAlmostEqual(sol.variables['xvars'][0].value, 1.0, places=places)
        self.assertAlmostEqual(sol.variables['xvars'][1].value, 4.743, places=places)
        self.assertAlmostEqual(sol.variables['xvars'][2].value, 3.82115, places=places)
        self.assertAlmostEqual(sol.variables['xvars'][3].value, 1.37941, places=places)

    def test_snopt(self):
        self.optimize('snopt')

    def test_slsqp(self):
        self.optimize('slsqp')

    def test_nlpqlp(self):
        self.optimize('nlpqlp')

    def test_fsqp(self):
        self.optimize('fsqp')

    def test_ipopt(self):
        self.optimize('ipopt')

    def test_conmin(self):
        opts = {'DELFUN' : 1e-9,
                'DABFUN' : 1e-9}
        self.optimize('conmin', optOptions=opts, places=2)

    def test_psqp(self):
        self.optimize('psqp')

if __name__ == "__main__":
    unittest.main()