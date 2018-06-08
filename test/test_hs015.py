"""Test solution of problem HS15 from the Hock & Schittkowski collection"""
from __future__ import print_function

import unittest

import numpy
from pyoptsparse import Optimization, OPT


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
        x = xdict['xvars']
        funcs = {}
        funcs['obj'] = [100*(x[1] - x[0]**2)**2 + (1-x[0])**2]
        conval = numpy.zeros(2, 'D')
        conval[0] = x[0]*x[1]
        conval[1] = x[0] + x[1]**2
        funcs['con'] = conval
        fail = False
        return funcs, fail

    def sens(self, xdict, funcs):
        x = xdict['xvars']
        funcsSens = {}
        funcsSens['obj'] = {'xvars': [2*100*(x[1]-x[0]**2)*(-2*x[0]) - 2*(1-x[0]),
                                      2*100*(x[1]-x[0]**2)]}
        funcsSens['con'] = {'xvars': [[x[1], x[0]],
                                      [1, 2*x[1]]]}
        fail = False

        return funcsSens, fail

    def optimize(self, optName, optOptions={}, storeHistory=False):
        # Optimization Object
        optProb = Optimization('HS15 Constraint Problem', self.objfunc)

        # Design Variables
        lower = [-5.0, -5.0]
        upper = [0.5,  5.0]
        value = [-2, 1.0]
        optProb.addVarGroup('xvars', 2, lower=lower, upper=upper, value=value)

        # Constraints
        lower = [1.0, 0.0]
        upper = [None, None]
        optProb.addConGroup('con', 2, lower=lower, upper=upper)

        # Objective
        optProb.addObj('obj')

        # Check optimization problem:
        # print(optProb)

        # Optimizer
        try:
            opt = OPT(optName, options=optOptions)
        except:
            raise unittest.SkipTest('Optimizer not available:', optName)

        # Solution
        if storeHistory:
            histFileName = '%s_hs015_Hist.hst' % (optName.lower())
        else:
            histFileName = None

        sol = opt(optProb, sens=self.sens, storeHistory=histFileName)

        # Check Solution
        self.assertAlmostEqual(sol.objectives['obj'].value, 306.5, places=5)

        self.assertAlmostEqual(sol.variables['xvars'][0].value, 0.5)
        self.assertAlmostEqual(sol.variables['xvars'][1].value, 2.0)

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
        self.optimize('conmin', optOptions=opts)

    def test_psqp(self):
        self.optimize('psqp')


if __name__ == "__main__":
    unittest.main()
