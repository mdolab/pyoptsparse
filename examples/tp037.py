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

def objfunc(xx):
    x = xx['xvars']
    fobj = -x[0]*x[1]*x[2]
    conval = [0]*2
    conval[0] = x[0] + 2.*x[1] + 2.*x[2] - 72.0
    conval[1] = -x[0] - 2.*x[1] - 2.*x[2]
    fcon = {'con':conval}
    fail = False

    return fobj, fcon, fail

# Optimization Object
optProb = Optimization('HS15 Constraint Problem', objfunc)

# Design Variables
optProb.addVarGroup('xvars',3, 'c',lower=[0,0,0], upper=[42,42,42], value=10)

# Constraints
optProb.addConGroup('con',2, lower=None, upper=0.0)

# Check optimization problem:
print optProb

# Optimizer
opt = OPT(options=optOptions)

# Solution
sol = opt(optProb, sens='CS')

# Check Solution
print sol

