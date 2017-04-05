# This example shows a bug where pySNOPT wouldn't optimize a model that has
# only equality constraints because it thought the problem was trivial. The
# problem is a simple paraboloid. The minimum should be at (7.166667,
# -7.833334), but with the bug, x and y stay at zero.

from __future__ import print_function
import numpy as np

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

def sens(xdict, funcs):
    """f(x,y) = (x-3)^2 + xy + (y+4)^2 - 3
    """
    x = xdict['x']
    y = xdict['y']
    funcsSens = {}

    funcsSens['obj', 'x'] = 2.0*x - 6.0 + y
    funcsSens['obj', 'y'] = 2.0*y + 8.0 + x
    #funcsSens['con', 'x'] = -1.0
    #funcsSens['con', 'y'] = 1.0

    fail = False
    return funcsSens, fail


con_jac = {}
con_jac['x'] = np.array(-1.0)
con_jac['y'] = np.array(1.0)


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
opt = SNOPT(optOptions = {'Major feasibility tolerance' : 1e-1})
sol = opt(optProb, sens=sens)

# Check Solution
print(sol)
print('Solution shoud be (x, y) = (7.166667, -7.833334)\n')