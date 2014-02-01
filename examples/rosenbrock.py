from __future__ import print_function
#!/usr/bin/env python
import time, sys
from pyoptsparse import Optimization
from pyoptsparse import SNOPT, IPOPT
import os, argparse
import numpy
import sys

parser = argparse.ArgumentParser()
parser.add_argument("--sens",help="sensitivity mode",type=str, default='FD')
parser.add_argument("--constrained",help="constrained or not",type=int,default=0)
parser.add_argument("--useDict",help="dictionary returns",type=int, default=0)
parser.add_argument("--testHist",help="test history",type=str,default="no")
parser.add_argument("--groups",help="use groups",type=int, default=0)
parser.add_argument("--sensMode",help="gradient mode",type=str, default='')

args = parser.parse_args()
sens = args.sens
constrained = args.constrained
useDict = args.useDict
testHist = args.testHist
groups = args.groups
sensMode = args.sensMode

def objfunc(xx):
    if groups:
        x = xx['x'] # Extract array
    else:
        x = xx

    fobj = 100*(x[1]-x[0]**2)**2+(1-x[0])**2

    if useDict:
        fcon = {}
        fcon['con'] = 0.1-(x[0]-1)**3 - (x[1]-1)
    else:
        if constrained:
            fcon = [0.1-(x[0]-1)**3 - (x[1]-1)]
        else:
            fcon = []

    fail = False

    return fobj, fcon, fail

def sensfunc(xx, fobj, fcon):
    if groups:
        x = xx['x'] # Extract array
    else:
        x = xx

    if useDict:
        gobj = {}
        gobj['xvars'] = [2*100*(x[1]-x[0]**2)*(-2*x[0]) - 2*(1-x[0]),
                     2*100*(x[1]-x[0]**2)]
    else:
        gobj = [2*100*(x[1]-x[0]**2)*(-2*x[0]) - 2*(1-x[0]),
                2*100*(x[1]-x[0]**2)]

    if useDict:
        gcon = {}
        gcon['con'] = {'xvars':[-3*(x[0]-1)**2, -1]}
    else:
        if constrained:
            gcon = [[-3*(x[0]-1)**2, -1]]
        else:
            gcon = []

    fail = False

    return gobj, gcon, fail

if sens == 'none':
    sens = None
if sens == 'user':
    sens = sensfunc

# Instantiate Optimization Problem
optProb = Optimization('Rosenbrock function', objfunc, useGroups=groups)
optProb.addVarGroup('x', 2, 'c', value=[-4,-4], lower=-5.12, upper=5.12,
                    scale=[1.0, 1.0], varSet='xvars')
if constrained:
    optProb.addCon('con',upper=0, scale=1.0)
optProb.addObj('f')

# Create optimizer
snopt = SNOPT()
ipopt = IPOPT()
if testHist == 'no':
    # Just run a normal run
    #solSnopt1 = snopt(optProb, sens=sens, sensMode=sensMode)#, storeHistory='opt_hist')
    solIpopt1 = ipopt(optProb, sens=sens, sensMode=sensMode)
    #print(solSnopt1.fStar)
    print(solIpopt1.fStar)
else:
    # First call just does 10 iterations
    snopt.setOption('Major iterations limit',10)
    solSnopt1 = snopt(optProb, sens=sens, sensMode='pgc', storeHistory='opt_hist')

    # Now we are allowed to do 50
    snopt.setOption('Major iterations limit',50)
    if testHist == 'hot':
        solSnopt2 = snopt(optProb, sens=sens, sensMode=sensMode,
                          hotStart='opt_hist', storeHistory='opt_hist')
    else:
        solSnopt2 = snopt(optProb, sens=sens, sensMode=sensMode,
                          coldStart='opt_hist', storeHistory='opt_hist')

    print(solSnopt2.fStar)
