from __future__ import print_function

import numpy
import argparse
from pyoptsparse import Optimization, OPT

parser = argparse.ArgumentParser()
parser.add_argument("--opt",help="optimizer",type=str, default='SLSQP')
args = parser.parse_args()
optOptions = {}

def objfunc(xx):
    x = xx['xvars']
    funcs = {}
    funcs['obj'] = -x[0]*x[1]*x[2]
    conval = [0]*2
    conval[0] = x[0] + 2.*x[1] + 2.*x[2] - 72.0
    conval[1] = -x[0] - 2.*x[1] - 2.*x[2]
    funcs['con'] = conval
    fail = False

    return funcs, fail

# Optimization Object
optProb = Optimization('TP037 Constraint Problem', objfunc)

# Design Variables
optProb.addVarGroup('xvars', 3, 'c',lower=[0,0,0], upper=[42,42,42], value=10)

# Constraints
optProb.addConGroup('con', 2, lower=None, upper=0.0)

# Objective
optProb.addObj('obj')

# Check optimization problem:
print(optProb)

# Optimizer
opt = OPT(args.opt, options=optOptions)

# Solution
sol = opt(optProb, sens='CS')

# Check Solution
print(sol)
