from __future__ import print_function
#!/usr/bin/env python
import time, sys
from pyoptsparse import Optimization, OPT
import os, argparse
import numpy
import sys
import pdb
# import cProfile

parser = argparse.ArgumentParser()
parser.add_argument("--sens",help="sensitivity mode",type=str, default='FD')
parser.add_argument("--constrained",help="constrained or not",type=int,default=0)
parser.add_argument("--testHist",help="test history",type=str,default="no")
parser.add_argument("--groups",help="use groups",type=int, default=0)
parser.add_argument("--sensMode",help="gradient mode",type=str, default='')
parser.add_argument("--opt",help="optimizer",type=str, default='SLSQP')
parser.add_argument('--pythonProf',help="profile in Python",type=int,default=0)
args = parser.parse_args()
sens = args.sens
constrained = args.constrained
testHist = args.testHist
groups = args.groups
sensMode = args.sensMode
pythonProf = args.pythonProf
optOptions = {}

def objfunc(xdict):
    x = xdict['xvars'] # Extract array
    funcs = {}
    funcs['obj'] = 100*(x[1]-x[0]**2)**2+(1-x[0])**2
    funcs['con'] = 0.1-(x[0]-1)**3 - (x[1]-1)
    fail = False

    return funcs, fail

def sensfunc(xdict, funcsDict):
    x = xdict['xvars'] # Extract array
    funcsSens = {}
    funcsSens['obj'] = {'xvars': [2*100*(x[1]-x[0]**2)*(-2*x[0]) - 2*(1-x[0]),
                                  2*100*(x[1]-x[0]**2)]}
    funcsSens['con'] = {'xvars':[-3*(x[0]-1)**2, -1]}
    fail = False
    return funcsSens, fail

# Matrix-free algorithm return functions
def objgrad(xdict):
    x = xdict['xvars']

    objsens = {}
    objsens['obj'] = {'xvars': numpy.array([2*100*(x[1]-x[0]**2)*(-2*x[0]) - 2*(1-x[0]),
                    2*100*(x[1]-x[0]**2)])}
    fail = False

    return objsens, fail

def jprod(xdict, pdict, sparse_only):
    x = xdict['xvars']
    p = pdict['xvars']
    qdict = {}

    if constrained:
        gcon = numpy.array([-3*(x[0]-1)**2, -1])
        qdict['con'] = numpy.array([numpy.dot(gcon,p)])
    else:
        q = numpy.zeros(0)
    fail = False

    return qdict, fail

def jtprod(xdict,qdict,sparse_only):
    x = xdict['xvars']

    pdict = {}
    if constrained:
        q = qdict['con']
        gcon = numpy.array([-3*(x[0]-1)**2, -1])
        pdict['xvars'] = q*gcon
    else:
        pdict['xvars'] = numpy.zeros(len(x))
    fail = False

    return pdict, fail

if sens == 'none':
    sens = None
if sens == 'user':
    sens = sensfunc
if sens == 'matrix-free':
    sens = [objgrad, jprod, jtprod]

# Instantiate Optimization Problem
optProb = Optimization('Rosenbrock function', objfunc)
optProb.addVarGroup('xvars', 2, 'c', value=[3, -3], lower=-5.12, upper=5.12,
                    scale=[1.0, 1.0])
optProb.finalizeDesignVariables()
if constrained:
    optProb.addCon('con',upper=0, scale=1.0)
optProb.addObj('obj')

# Create optimizer
opt = OPT(args.opt, options=optOptions)
if testHist == 'no':
    # Just run a normal run
    sol = opt(optProb, sens=sens, sensMode=sensMode)
    # print(sol.fStar)
    print(sol)
else:
    # First call just does 10 iterations
    if args.opt.lower() == 'snopt':
        opt.setOption('Major iterations limit',10)
        solSnopt1 = opt(optProb, sens=sens, sensMode='pgc', storeHistory='opt_hist')

        # Now we are allowed to do 50
        opt.setOption('Major iterations limit',50)
        if testHist == 'hot':
            solSnopt2 = opt(optProb, sens=sens, sensMode=sensMode,
                              hotStart='opt_hist', storeHistory='opt_hist')
        else:
            solSnopt2 = opt(optProb, sens=sens, sensMode=sensMode,
                              coldStart='opt_hist', storeHistory='opt_hist')

        print(solSnopt2.fStar)

    else:
        # Quick test that history prints ok
        # if pythonProf:
        #     cProfile.run('sol = opt(optProb, sens=sens, sensMode=sensMode, storeHistory=\'opt_hist.bin\')',args.opt.lower()+'_profile.profile')
        # else:
        sol = opt(optProb, sens=sens, sensMode=sensMode, storeHistory='opt_hist.bin')
