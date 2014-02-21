#!/usr/bin/env python
import numpy
import argparse
from pyoptsparse import Optimization

parser = argparse.ArgumentParser()
parser.add_argument("--opt",help="optimizer",type=str, default='SNOPT')
args = parser.parse_args()

optOptions = {}
if args.opt.lower() == 'ipopt':
    from pyoptsparse import IPOPT as OPT
elif args.opt.lower() == 'snopt':
    from pyoptsparse import SNOPT as OPT
elif args.opt.lower() == 'slsqp':
    from pyoptsparse import SLSQP as OPT
elif args.opt.lower() == 'conmin':
    from pyoptsparse import CONMIN as OPT
elif args.opt.lower() == 'fsqp':
    from pyoptsparse import FSQP as OPT
elif args.opt.lower() == 'nlpql':
    from pyoptsparse import NLPQL as OPT

def objfunc(xdict):
    x = xdict['xvars'] # Extract array
    fobj = 100*(x[1]-x[0]**2)**2+(1-x[0])**2
    fcon = {}
    fail = False

    return fobj, fcon, fail

def sens(xdict, fobj, fcon):
    x = xdict['xvars'] # Extract array
    gobj = {}
    gobj['xvars'] = [2*100*(x[1]-x[0]**2)*(-2*x[0]) - 2*(1-x[0]),
                     2*100*(x[1]-x[0]**2)]
    gcon = {}
    fail = False

    return gobj, gcon, fail


# Optimization Object
optProb = Optimization('TP109 Constraint Problem',objfunc)

# Design Variables
optProb.addVarGroup('xvars', 2, value=0)

# Constraints -- None

# Check optimization problem:
print optProb
optProb.printSparsity()

# Optimizer
opt = OPT(options=optOptions)

# Solution
sol = opt(optProb, sens=sens, storeHistory='opt_hist')

# Check Solution
print sol
