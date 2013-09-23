#!/usr/bin/env python

import os,sys

from pyOpt_constraint import Constraint
from pyOpt_optimizer import Optimizer
from pyOpt_optimization import Optimization

__all__ = ['Constraint','Optimization']

dir = os.path.dirname(os.path.realpath(__file__))
for f in os.listdir(dir):
    if f.startswith('py') and os.path.isdir(os.path.join(dir,f)):
        try:
            exec 'from %s import %s' %(f,f.strip('py'))
            __all__.extend(sys.modules['pyOpt.'+f].__all__)
        except:
            continue
        #end
    #end
#end
