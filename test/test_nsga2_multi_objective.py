""" Test NSGA2."""
from __future__ import print_function

import unittest
from numpy.testing import assert_allclose
from pyoptsparse import Optimization, NSGA2


def objfunc(xdict):
    x = xdict['x']
    y = xdict['y']

    funcs = {}
    funcs['obj1'] = (x - 0.0)**2 + (y - 0.0)**2
    funcs['obj2'] = (x - 1.0)**2 + (y - 1.0)**2

    fail = False

    return funcs, fail


class TestNSGA2(unittest.TestCase):

    def test_opt(self):
        # Instantiate Optimization Problem
        optProb = Optimization('Rosenbrock function', objfunc)
        optProb.addVar('x', 'c', value=0, lower=-600, upper=600)
        optProb.addVar('y', 'c', value=0, lower=-600, upper=600)

        optProb.addObj('obj1')
        optProb.addObj('obj2')

        #300 generations will find x=(0,0), 200 or less will find x=(1,1)
        options = {
            'maxGen':200
        }

        # Optimizer
        try:
            opt = NSGA2(options=options)
        except:
            raise unittest.SkipTest('Optimizer not available:', 'NSGA2')

        sol = opt(optProb)

        # Check Solution
        tol = 1E-2
        assert_allclose(sol.variables['x'][0].value, 1.0, atol=tol, rtol=tol)

        assert_allclose(sol.variables['y'][0].value, 1.0, atol=tol, rtol=tol)


if __name__ == "__main__":
    unittest.main()