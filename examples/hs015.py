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

import numpy
import argparse
from pyoptsparse import Optimization, OPT
parser = argparse.ArgumentParser()
parser.add_argument("--opt", help="optimizer", type=str, default='SNOPT')
parser.add_argument("--storeHistory",help="option to store history",type=int, default=0)
args = parser.parse_args()
optOptions = {}

def objfunc(xdict):
    x = xdict['xvars']
    fobj = 100*(x[1] - x[0]**2)**2 + (1-x[0])**2
    conval = numpy.zeros(2, 'D')
    conval[0] = x[0]*x[1]
    conval[1] = x[0] + x[1]**2
    fcon = {'con':conval}

    fail = False

    return fobj, fcon, fail

def sens(xdict, fobj, fcon):
    x = xdict['xvars']
    gobj = {}
    gobj['xvars'] = [2*100*(x[1]-x[0]**2)*(-2*x[0]) - 2*(1-x[0]),
                     2*100*(x[1]-x[0]**2)]
    gcon = {}
    gcon['con'] = {'xvars':[[x[1], x[0]],
                            [1   , 2*x[1]]]}
    fail = False
    
    return gobj, gcon, fail

# Optimization Object
optProb = Optimization('HS15 Constraint Problem', objfunc)

# Design Variables
lower = [-5, -5]
upper = [0.5,  5]
value = [-2, 1]
optProb.addVarGroup('xvars', 2, lower=lower, upper=upper, value=value)

# Constraints
lower = [1, 0]
upper = [None, None]
optProb.addConGroup('con', 2, lower=lower, upper=upper)

# Check optimization problem:
print optProb

# Optimizer
opt = OPT(args.opt, options=optOptions)

# Solution
if args.storeHistory:
    histFileName = '%s_hs015_Hist.hst'%(args.opt.lower())
else:
    histFileName = None
# end
sol = opt(optProb, sens=sens,storeHistory=histFileName)

# Check Solution
print sol
