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
    x = xdict['xvars']
    fobj = x[0] * x[3] * (x[0] + x[1] + x[2]) + x[2]
    fcon = {}
    fcon['con'] = [x[0] * x[1] * x[2] * x[3], 
                   x[0]*x[0] + x[1]*x[1] + x[2]*x[2] + x[3]*x[3]]
    fail = False
    return fobj, fcon, False

def sens(xdict, fobj, fcon):
    x = xdict['xvars']
    gobj = {};
    gobj = numpy.array([x[0] * x[3] + x[3] * (x[0] + x[1] + x[2]) , 
                        x[0] * x[3],
                        x[0] * x[3] + 1.0,
                        x[0] * (x[0] + x[1] + x[2])
                        ])
    gcon = {}
    jac = [[x[1]*x[2]*x[3], x[0]*x[2]*x[3], x[0]*x[1]*x[3], x[0]*x[1]*x[2]],
           [2.0*x[0], 2.0*x[1], 2.0*x[2], 2.0*x[3]]]
    gcon['con'] = {'xvars': jac}
    fail = False
    return gobj, gcon, fail

# Optimization Object
optProb = Optimization('HS071 Constraint Problem', objfunc)

# Design Variables
x0 = [1.0, 5.0, 5.0, 1.0]
optProb.addVarGroup('xvars', 4, lower=1, upper=5, value=x0)

# Constraints
optProb.addConGroup('con', 2, lower=[25,40], upper=[1e19, 40])

# Check optimization problem:
print optProb

# Optimizer
opt = OPT(options=optOptions)

# Solution
sol = opt(optProb, sens=sens)

# Check Solution
print sol

