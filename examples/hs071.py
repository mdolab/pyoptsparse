import numpy
import argparse
from pyoptsparse import Optimization, OPT

parser = argparse.ArgumentParser()
parser.add_argument("--opt",help="optimizer",type=str, default='SNOPT')
args = parser.parse_args()
optOptions = {}

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

# Optimization Object
optProb = Optimization('HS071 Constraint Problem', objfunc)

# Design Variables
x0 = [1.0, 5.0, 5.0, 1.0]
optProb.addVarGroup('xvars', 4, lower=1, upper=5, value=x0)

# Constraints
optProb.addConGroup('con', 2, lower=[25,40], upper=[1e19, 40])

# Objective
optProb.addObj('obj')

# Check optimization problem:
print optProb

# Optimizer
opt = OPT(args.opt, options=optOptions)

# Solution
sol = opt(optProb, sens=sens)

# Check Solution
print sol

