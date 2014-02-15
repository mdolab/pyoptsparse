#!/usr/bin/env python
"""
Solves Schittkowski's TP109 constraint problem.

    min 	3.0*x1+1.*10**(-6)*x1**3+0.522074*10**(-6)*x2**3+2.0*x2 
    s.t.    -(x4-x3+0.550) <= 0
            -(x3-x4+0.550) <= 0
            -(2.25*10**(+6)-x1**2-x8**2) <= 0
            -(2.25*10**(+6)-x2**2-x9**2) <= 0
            -(x5*x6*sin(-x3-0.250) + +x5*x7*sin(-x4-0.250)+2.0*x5**2*b)*ra+400.0-x1 = 0
            -((x5*x6*sin(x3-0.250)+x6*x7*sin(x3-x4-0.250)+2.0*x6**2*b)*ra+400.0-x2) = 0
            -((x5*x7*sin(x4-0.250)+x6*x7*sin(x4-x3-0.250)+2.0*x7**2*b)*ra+881.7790) = 0
            -(x8+(x5*x6*cos(-x3-0.250)+x5*x7*cos(-x4-0.250)-2.0*x5**2*c)*ra+0.7533*10**(-3)*x5**2-200.00) = 0
            -(x9+(x5*x6*cos(x3-0.250)+x7*x6*cos(x3-x4-0.250)-2.0*x6**2*c)*ra+0.7533*10**(-3)*x(6)**2-200.00) = 0
            -((x5*x7*cos(x4-0.250)+x6*x7*cos(x4-x3-0.250)-2.0*x7**2*c)*ra+0.7533*10**(-3)*x7**2-22.9380) = 0
            0 <= xi, i = 1,2
            -0.55 <= xi <= 0.55, i = 3,4
            196.0 <= xi <= 252.0, i = 5,6,7
            -400.0 <= xi <= 800.0, i = 8,9
    where   a = 50.176
            b = sin(0.25)
            c = cos(0.25)
            ra = 1.0/50.176

            
    f*1 = 0.536206927538e+04
    x*1 = [0.674888100445e+03, 0.113417039470e+04, 0.133569060261e+00, -0.371152592466e+00, 0.252e+03, 0.252e+03, 0.201464535316e+03, 0.426660777226e+03, 0.368494083867e+03]
"""
import numpy, sys
from numpy import sin, cos
from pyoptsparse import Optimization
if 'IPOPT' in sys.argv:
    from pyoptsparse import IPOPT as OPT
    optOptions={'max_iter':150,
                'tol': 1e-6,
                'derivative_test':'first-order',
                'derivative_test_print_all':'no',
                'derivative_test_tol':1e-4,
                'output_file':'testoutIPOPT.out'}
elif 'SLSQP' in sys.argv:
    from pyoptsparse import SLSQP as OPT
    optOptions = {}
elif 'CONMIN' in sys.argv:
    from pyoptsparse import CONMIN as OPT
    optOptions = {}
elif 'FSQP' in sys.argv:
    from pyoptsparse import FSQP as OPT
    optOptions = {}
elif 'NLPQL' in sys.argv:
    from pyoptsparse import NLPQL as OPT
    optOptions = {}
else:
    from pyoptsparse import SNOPT as OPT
    optOptions={}

USE_LINEAR = False
if 'linear' in sys.argv:
    USE_LINEAR=True

def objfunc(xx):
    x = xx['xvars']
    
    a=50.1760
    b=sin(0.250)    
    c=cos(0.250)

    fobj = 3.0*x[0] + (1e-6)*x[0]**3 + 0.522074e-6*x[1]**3 + 2*x[1]
    fcon = numpy.zeros(10,'D')
    fcon[0] = 2250000 - x[0]**2 - x[7]**2

    fcon[1] = 2250000 - x[1]**2 - x[8]**2

    fcon[2] = x[4]*x[5]*sin(-x[2]-0.25) + x[4]*x[6]*sin(-x[3] - 0.25) + 2*b*x[4]**2 - a*x[0] + 400*a
    
    fcon[3] = x[4]*x[5]*sin(x[2] - 0.25) + x[5]*x[6]*sin(x[2] - x[3] - 0.25) + 2*b*x[5]**2 - a*x[1] + 400*a

    fcon[4] = x[4]*x[6]*sin(x[3] - 0.25) + x[5]*x[6]*sin(x[3] - x[2] - 0.25) + 2*b*x[6]**2 + 881.779*a

    fcon[5] = a*x[7] + x[4]*x[5]*cos(-x[2] - 0.25) + x[4]*x[6]*cos(-x[3] - 0.25) - 200*a - 2*c*x[4]**2 + 0.7533e-3*a*x[4]**2

    fcon[6] = a*x[8] + x[4]*x[5]*cos(x[2]-0.25) + x[5]*x[6]*cos(x[2] - x[3] - 0.25) - 2*c*x[5]**2 + 0.7533e-3*a*x[5]**2 - 200*a

    fcon[7] = x[4]*x[6]*cos(x[3] - 0.25) + x[5]*x[6]*cos(x[3] - x[2] - 0.25) - 2*c*x[6]**2 - 22.938*a + 0.7533e-3*a*x[6]**2

    fcon[8] = x[3] - x[2] + 0.55
    fcon[9] = x[2] - x[3] + 0.55

    if USE_LINEAR:
        fcon = {'con':fcon[0:8]}
    else:
        fcon = {'con':fcon[0:10]}
    fail = False
  
    return fobj, fcon, fail
# 
# ============================================================================= 
optProb = Optimization('TP109 Constraint Problem',objfunc)
lower = [0.0, 0.0, -0.55, -0.55, 196, 196, 196, -400, -400]
upper = [None, None, 0.55, 0.55,  252, 252, 252, 800, 800]
value = [0,0,0,0,0,0,0,0,0]

#value_opt = [0.7e+03, 0.1e+04, 0.1, -0.3, 200, 200, 200, 0.4e+03, 0.4e+03]
#value=value_opt
optProb.addVarGroup('xvars', 9, lower=lower, upper=upper, value=value)

lower = [0, 0      , 0, 0, 0, 0, 0, 0]
upper = [None, None, 0, 0, 0, 0, 0, 0]
if not USE_LINEAR:
    lower.extend([0,0])
    upper.extend([None, None])

optProb.finalizeDesignVariables()
# We separate out the 8 non-linear constraints

optProb.addConGroup('con', len(lower), lower=lower, upper=upper)
# And the 2 linear constriants
if USE_LINEAR:
    # jac = numpy.zeros((2, 9))
    # jac[0, 3] = 1.0; jac[0, 2] = -1.0
    # jac[1, 2] = 1.0; jac[1, 3] = -1.0
    # optProb.addConGroup('lin_con',2, lower=[-.55,-.55], upper=[None, None],
    #                     wrt=['xvars'], jac ={'xvars':jac}, linear=True)

# OR do it with a single two-sided constraint
    jac = numpy.zeros((1, 9))
    jac[0, 3] = 1.0; jac[0, 2] = -1.0
    optProb.addConGroup('lin_con',1, lower=[-.55], upper=[0.55],
                        wrt=['xvars'], jac ={'xvars':jac}, linear=True)



print optProb
optProb.addObj('f')
optProb.printSparsity()
opt = OPT(options=optOptions)
sol = opt(optProb, sens='CS')
#print sol
