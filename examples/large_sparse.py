# This example is taken from the OpenOpt Examples website. 
# http://trac.openopt.org/openopt/browser/PythonPackages/FuncDesigner/FuncDesigner/examples/nlpSparse.py
from __future__ import print_function

import numpy
from numpy import arange
from scipy import sparse
import argparse
from pyoptsparse import Optimization, OPT

parser = argparse.ArgumentParser()
parser.add_argument("--opt",help="optimizer",type=str, default='SNOPT')
parser.add_argument('--optOptions', type=str, help='Options for the optimizer', default="{}")
args = parser.parse_args()
exec('optOptions=%s'% args.optOptions)

def objfunc(xdict):
    x = xdict['x']; y = xdict['y']; z = xdict['z']
    funcs = {}
    funcs['obj'] = x**2 + 2*numpy.sum(y**2) + 3*numpy.sum(z)
    funcs['con1'] = x + 1e-3*abs(x)**2.05
    funcs['con2'] = x**4 + numpy.sum(y) + numpy.sum(z**2)
    funcs['con3'] = x + numpy.sum(z)
    
    return funcs, False

def sens(xdict, funcs):
    x = xdict['x'];  y = xdict['y']; z = xdict['z']
    funcsSens = {}
    funcsSens = {('obj', 'x') : [2*x],
                 ('obj', 'y') : 4*y,
                 ('obj', 'z') : 3*numpy.ones(2*N),
                 ('con1', 'x') : 2.05*x*(x*x)**0.025,
                 ('con2', 'x') : 4*x**3,
                 ('con2', 'y') : numpy.ones(N),
                 ('con2', 'z') : 2*z,
                 ('con3', 'x') : 1.0,
                 ('con3', 'z') : numpy.ones(2*N)}
    
    return funcsSens, False

# Optimization Object
optProb = Optimization('large and sparse', objfunc)

N = 50000
# Design Variables
optProb.addVar('x', lower=-100, upper=150, value=0)
optProb.addVarGroup('y', N, lower=-10-arange(N), upper=arange(N), value=0)
optProb.addVarGroup('z', 2*N, upper=arange(2*N), lower=-100-arange(2*N), value=0)
# Constraints
optProb.addCon('con1', upper=100, wrt=['x'])
optProb.addCon('con2', upper=100)
optProb.addCon('con3', lower=4, wrt=['x','z'])
optProb.addConGroup('lincon', N, lower=2-3*arange(N), linear=True,
                    wrt=['x','y'], 
                    jac={'x': numpy.ones((N,1)),
                         'y':sparse.spdiags(numpy.ones(N), 0, N, N)})
optProb.addObj('obj')
# Optimizer
opt = OPT(args.opt, options=optOptions)
optProb.printSparsity()
# Solution
sol = opt(optProb, sens=sens)
print(sol)
