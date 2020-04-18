# This example shows a bug where pySNOPT wouldn't optimize a model that has
# only equality constraints because it thought the problem was trivial. The
# problem is a simple paraboloid. The minimum should be at (7.166667,
# -7.833334), but with the bug, x and y stay at zero.

from __future__ import print_function
import unittest

import numpy as np
from numpy.testing import assert_almost_equal
from pyoptsparse import Optimization, SNOPT


def objfunc(xdict):
    """ Evaluates the equation f(x,y) = (x-3)^2 + xy + (y+4)^2 - 3 """
    x = xdict['x']
    y = xdict['y']
    funcs = {}

    funcs['obj'] =  (x-3.0)**2 + x*y + (y+4.0)**2 - 3.0
    conval = -x + y
    funcs['con'] = conval

    fail = False
    return funcs, fail

def objfunc_no_con(xdict):
    """ Evaluates the equation f(x,y) = (x-3)^2 + xy + (y+4)^2 - 3 """
    x = xdict['x']
    y = xdict['y']
    funcs = {}

    funcs['obj'] =  (x-3.0)**2 + x*y + (y+4.0)**2 - 3.0

    fail = False
    return funcs, fail

def objfunc_2con(xdict):
    """ Evaluates the equation f(x,y) = (x-3)^2 + xy + (y+4)^2 - 3 """
    x = xdict['x']
    y = xdict['y']
    funcs = {}

    funcs['obj'] =  (x-3.0)**2 + x*y + (y+4.0)**2 - 3.0
    conval = -x + y
    funcs['con'] = conval *np.ones(2)
    funcs['con2'] = (conval + 1) *np.ones(3)

    fail = False
    return funcs, fail

def sens(xdict, funcs):
    """f(x,y) = (x-3)^2 + xy + (y+4)^2 - 3
    """
    x = xdict['x']
    y = xdict['y']
    funcsSens = {}

    funcsSens['obj', 'x'] = 2.0*x - 6.0 + y
    funcsSens['obj', 'y'] = 2.0*y + 8.0 + x

    fail = False
    return funcsSens, fail


con_jac = {}
con_jac['x'] = np.array(-1.0)
con_jac['y'] = np.array(1.0)


class TestSNOPTBug(unittest.TestCase):

    def test_opt(self):
        # Optimization Object
        optProb = Optimization('Paraboloid', objfunc)

        # Design Variables
        optProb.addVarGroup('x', 1, type='c', lower=-50.0, upper=50.0, value=0.0)
        optProb.addVarGroup('y', 1, type='c', lower=-50.0, upper=50.0, value=0.0)
        optProb.finalizeDesignVariables()

        # Objective
        optProb.addObj('obj')

        # Equality Constraint
        optProb.addConGroup('con', 1, lower=-15.0, upper=-15.0, wrt=['x', 'y'], linear=True, jac=con_jac)

        # Check optimization problem:
        print(optProb)


        # Optimizer
        try:
            opt = SNOPT(optOptions = {'Major feasibility tolerance' : 1e-1})
        except:
            raise unittest.SkipTest('Optimizer not available: SNOPT')

        sol = opt(optProb, sens=sens)

        # Check Solution 7.166667, -7.833334
        assert_almost_equal(sol.variables['x'][0].value, 7.166667, decimal=6)
        assert_almost_equal(sol.variables['y'][0].value, -7.833333, decimal=6)

    def test_opt_bug1(self):
        # Due to a new feature, there is a TypeError when you optimize a model without a constraint.
        optProb = Optimization('Paraboloid', objfunc_no_con)

        # Design Variables
        optProb.addVarGroup('x', 1, type='c', lower=-50.0, upper=50.0, value=0.0)
        optProb.addVarGroup('y', 1, type='c', lower=-50.0, upper=50.0, value=0.0)
        optProb.finalizeDesignVariables()

        # Objective
        optProb.addObj('obj')

        # Optimizer
        try:
            opt = SNOPT(optOptions = {'Major feasibility tolerance' : 1e-1})
        except:
            raise unittest.SkipTest('Optimizer not available: SNOPT')

        sol = opt(optProb, sens=sens)

    def test_opt_bug_print_2con(self):
        # Optimization Object
        optProb = Optimization('Paraboloid', objfunc_2con)

        # Design Variables
        optProb.addVarGroup('x', 1, type='c', lower=-50.0, upper=50.0, value=0.0)
        optProb.addVarGroup('y', 1, type='c', lower=-50.0, upper=50.0, value=0.0)
        optProb.finalizeDesignVariables()

        # Objective
        optProb.addObj('obj')

        con_jac2 = {}
        con_jac2['x'] = -np.ones((2, 1))
        con_jac2['y'] = np.ones((2, 1))

        con_jac3 = {}
        con_jac3['x'] = -np.ones((3, 1))
        con_jac3['y'] = np.ones((3, 1))

        # Equality Constraint
        optProb.addConGroup('con', 2, lower=-15.0, upper=-15.0, wrt=['x', 'y'], linear=True, jac=con_jac2)
        optProb.addConGroup('con2', 3, lower=-15.0, upper=-15.0, wrt=['x', 'y'], linear=True, jac=con_jac3)

        # Check optimization problem:
        print(optProb)

        # Optimizer
        try:
            opt = SNOPT(optOptions = {'Major feasibility tolerance' : 1e-1})
        except:
            raise unittest.SkipTest('Optimizer not available: SNOPT')

        sol = opt(optProb, sens=sens)

        print(sol)

if __name__ == "__main__":
    unittest.main()