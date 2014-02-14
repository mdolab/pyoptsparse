from __future__ import absolute_import
#!/usr/bin/env python
#from . import pyOptSparse

# By default we will import the constraint and optProb classes since
# these are required for all optimizers. All other optimizers will be
# have to be imported manually 
__all__ = []

# Import all the core modules
from .pyOptSparse.pyOpt_history import History
from .pyOptSparse.pyOpt_variable import Variable
from .pyOptSparse.pyOpt_constraint import Constraint
from .pyOptSparse.pyOpt_objective import Objective
from .pyOptSparse.pyOpt_optimization import Optimization
from .pyOptSparse.pyOpt_optimizer import Optimizer
from .pyOptSparse.pyOpt_gradient import Gradient
from .pyOptSparse.pyOpt_solution import Solution
# All defines for * imports
__all__.append('Constraint')
__all__.append('Optimization')

# Now try to import all the individual optimizers
try:
    from .pyOptSparse.pySNOPT.pySNOPT import SNOPT
except ImportError:
    pass

try:
    from .pyOptSparse.pyIPOPT.pyIPOPT import IPOPT
except ImportError:
    pass

try:
    from .pyOptSparse.pySLSQP.pySLSQP import SLSQP
except ImportError:
    pass

try:
    from .pyOptSparse.pyCONMIN.pyCONMIN import CONMIN
except ImportError:
    pass


# try:
#     from .pyOptSparse.pyNSGA2.pyNSGA2 import NSGA2
# except ImportError:
#     pass

