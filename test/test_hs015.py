"""Test solution of problem HS15 from the Hock & Schittkowski collection"""
from __future__ import print_function

import unittest

import numpy as np
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
        self.nf += 1
        x = xdict['xvars']
        funcs = {}
        funcs['obj'] = [100*(x[1] - x[0]**2)**2 + (1-x[0])**2]
        conval = np.zeros(2, 'D')
        conval[0] = x[0]*x[1]
        conval[1] = x[0] + x[1]**2
        funcs['con'] = conval
        fail = False
        return funcs, fail

    def sens(self, xdict, funcs):
        self.ng += 1
        x = xdict['xvars']
        funcsSens = {}
        funcsSens['obj'] = {'xvars': [2*100*(x[1]-x[0]**2)*(-2*x[0]) - 2*(1-x[0]),
                                      2*100*(x[1]-x[0]**2)]}
        funcsSens['con'] = {'xvars': [[x[1], x[0]],
                                      [1, 2*x[1]]]}
        fail = False
        return funcsSens, fail

    def optimize(self, optName, optOptions={}, storeHistory=False,places=5, hotStart=None):
        self.nf = 0 # number of function evaluations
        self.ng = 0 # number of gradient evaluations
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
            self.histFileName = '%s_hs015_Hist.hst' % (optName.lower())
        else:
            self.histFileName = None

        sol = opt(optProb, sens=self.sens, storeHistory=self.histFileName, hotStart=hotStart)

        # Test printing solution to screen
        print(sol)

        # Check Solution
        fobj = sol.objectives['obj'].value
        diff = np.min(np.abs([fobj - 306.5, fobj - 360.379767]))
        self.assertAlmostEqual(diff, 0.0, places=places)

        xstar1 = (0.5, 2.0)
        xstar2 = (-0.79212322, -1.26242985)
        x1 = sol.variables['xvars'][0].value
        x2 = sol.variables['xvars'][1].value

        dv = sol.getDVs()
        self.assertAlmostEqual(x1, dv['xvars'][0], places=10)
        self.assertAlmostEqual(x2, dv['xvars'][1], places=10)

        diff = np.min(np.abs([xstar1[0] - x1, xstar2[0] - x1]))
        self.assertAlmostEqual(diff, 0.0, places=places)

        diff = np.min(np.abs([xstar1[1] - x2, xstar2[1] - x2]))
        self.assertAlmostEqual(diff, 0.0, places=places)

    def test_snopt(self):
        store_vars = ['step','merit','feasibility','optimality','penalty','Hessian','condZHZ','slack','lambda']
        optOptions = {
            'Save major iteration variables': store_vars
        }
        self.optimize('snopt',optOptions=optOptions,storeHistory=True)
        from sqlitedict import SqliteDict
        hist = SqliteDict(self.histFileName)
        self.assertIn('isMajor',hist['0'].keys())
        self.assertEqual(7,hist['19']['nMajor'])
        for var in store_vars:
            self.assertIn(var,hist['19'].keys())
        # re-optimize with hotstart
        self.optimize('snopt',storeHistory=False,hotStart=self.histFileName)
        # now we should do the same optimization without calling them
        self.assertEqual(self.nf,0)
        self.assertEqual(self.ng,0)

    def test_slsqp(self):
        self.optimize('slsqp', storeHistory=True)
        self.assertGreater(self.nf,0)
        self.assertGreater(self.ng,0)
        # re-optimize with hotstart
        self.optimize('slsqp',storeHistory=False,hotStart=self.histFileName)
        # now we should do the same optimization without calling them
        self.assertEqual(self.nf,0)
        self.assertEqual(self.ng,0)


    def test_nlpqlp(self):
        self.optimize('nlpqlp', storeHistory=True)
        self.assertGreater(self.nf,0)
        self.assertGreater(self.ng,0)
        # re-optimize with hotstart
        self.optimize('nlpqlp',storeHistory=False,hotStart=self.histFileName)
        # now we should do the same optimization without calling them
        self.assertEqual(self.nf,0)
        self.assertEqual(self.ng,0)

    def test_ipopt(self):
        self.optimize('ipopt', places=4, storeHistory=True)
        self.assertGreater(self.nf,0)
        self.assertGreater(self.ng,0)
        # re-optimize with hotstart
        self.optimize('ipopt',storeHistory=False,hotStart=self.histFileName, places=4)
        # now we should do the same optimization without calling them
        self.assertEqual(self.nf,0)
        self.assertEqual(self.ng,0)

    def test_paropt(self):
        self.optimize('paropt', storeHistory=True)
        self.assertGreater(self.nf,0)
        self.assertGreater(self.ng,0)
        # re-optimize with hotstart
        self.optimize('paropt',storeHistory=False,hotStart=self.histFileName)
        # now we should do the same optimization without calling them
        self.assertEqual(self.nf,0)
        self.assertEqual(self.ng,0)

    def test_conmin(self):
        opts = {'DELFUN' : 1e-9,
                'DABFUN' : 1e-9}
        self.optimize('conmin', optOptions=opts, storeHistory=True)
        # re-optimize with hotstart
        self.optimize('conmin',optOptions=opts,storeHistory=False,hotStart=self.histFileName)
        # now we should do the same optimization without calling them
        self.assertEqual(self.nf,0)
        self.assertEqual(self.ng,0)

    def test_psqp(self):
        self.optimize('psqp', storeHistory=True)
        self.assertGreater(self.nf,0)
        self.assertGreater(self.ng,0)
        # re-optimize with hotstart
        self.optimize('psqp',storeHistory=False,hotStart=self.histFileName)
        # now we should do the same optimization without calling them
        self.assertEqual(self.nf,0)
        self.assertEqual(self.ng,0)

if __name__ == "__main__":
    unittest.main()
