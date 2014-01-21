#!/usr/bin/env python
"""
pyOpt_solution 

This class is used to describe the solution of an optimization
problem. This class is inherits from Optimization which enables a
solution to be used as an input to a subsequent optimization problem.

Copyright (c) 2013-2014 by Dr. Gaetan Kenway
All rights reserved.

Developers:
-----------
- Dr. Gaetan K.W. Kenway (GKK)

History
-------
    v. 1.0  - Initial Class Creation (GKK, 2013)
"""

# =============================================================================
# Standard Python modules
# =============================================================================
import copy
from pyOpt_optimization import Optimization
class Solution(Optimization):
    """
    Solution Class Initialization
    
    Parameters
    ----------
    optProb : Optimization problem class
        Optimization problem used to create solution

    optTime : float
        Time required for the optimziation
        
    optEvals : int
        The number of function evalution for the solution

    optInform : int
        The inform code from the optimization.
        """
    
    def __init__(self, optProb, optTime, optEvals, optInform):
        
        Optimization.__init__(self, optProb.name, optProb.objFun, optProb.useGroups)

        # Copy over the variables, constraints, and objectives
        self.variables = copy.deepcopy(optProb.variables)
        self.constraints = copy.deepcopy(optProb.constraints)
        self.objectives = copy.deepcopy(optProb.objectives)
        self.varSetNames = copy.deepcopy(optProb.varSetNames)
        self.varGroupNames = copy.deepcopy(optProb.varGroupNames)
        self.conGroupNames = copy.deepcopy(optProb.conGroupNames)
        self.optTime = optTime
        self.optEvals = optEvals
        self.optInform = optInform

    def __str__(self):
        """
        Print Structured Solution
        
        Documentation last updated:  April. 30, 2008 - Peter W. Jansen
        """
        
        text0 = Optimization.__str__(self)
        text1 = ''
        lines = text0.split('\n')
        lines[1] = lines[1][len('Optimization Problem -- '):]
        for i in xrange(5):
            text1 += lines[i] + '\n'
        #end
        # if self.display_opt:
        #     text1 += '\n	Options:\n '
        #     opt_keys = self.options_set.keys()
        #     opt_keys.sort()
        #     for key in opt_keys:
        #         ns = 25-len(key)
        #         text1 += '		'+ key +':' + str(self.options_set[key][1]).rjust(ns,'.') + '\n'
        #     #end
        # #end
        text1 += '\n    Solution: \n'
        text1 += ('-'*80) + '\n'
        text1 += '    Total Time: %25.4f\n' %(self.optTime)
        text1 += '    Total Function Evaluations: %9.0i\n' %(self.optEvals)
        # for key in self.parameters.keys():
        #     if (isinstance(self.parameters[key],(dict,list,tuple))) and (len(self.parameters[key]) == 0):
        #         continue
        #     elif (isinstance(self.parameters[key],numpy.ndarray)) and (0 in (self.parameters[key]).shape):
        #         continue
        #     else:
        #         text1 += '    '+ key +': ' + str(self.parameters[key]).rjust(9) + '\n'
        #     #end
        # #end

        for i in xrange(5,len(lines)):
            text1 += lines[i] + '\n'
        #end
        text1 += ('-'*80) + '\n'

        return text1
