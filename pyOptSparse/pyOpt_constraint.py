#!/usr/bin/env python
"""
pyOpt_constraint

Holds the representation of a pyOptSparse constraint group

Copyright (c) 2008-2013 by pyOpt Developers
All rights reserved.

Developers:
-----------
- Dr. Gaetan K.W. Kenway (GKK)

History
-------
    v. 1.0  - Initial Class Creation (GKK, 2013)
"""
# =============================================================================
# External Python modules
# =============================================================================
import numpy
#from pyOpt_optimization import INFINITY
INFINITY=1e20
eps = numpy.finfo(1.0).eps
# =============================================================================
# Constraint Class
# =============================================================================
class Constraint(object):
    """
    Constraint Class Initialization
    """
    def __init__(self, name, nCon, linear, wrt, jac, partialReturnOk,
                 lower, upper, scale):

        self.name = name
        self.ncon = nCon
        self.linear = linear
        self.wrt = wrt
        self.jac = jac
        self.partialReturnOk = partialReturnOk

        # Indices of this constraint. rs='row start', re='row end'. 
        self.rs = None
        self.re = None
        self.scale = scale # Unused, just store for reference

        # Before we can do the processing below we need to have lower
        # and upper arguments expanded:

        if lower is None:
            lower = [None for i in range(self.ncon)]
        elif numpy.isscalar(lower):
            lower = lower*numpy.ones(self.ncon)
        elif len(lower) == self.ncon:
            pass  # Some iterable object
        else:
            raise Error('The \'lower\' argument to addCon or addConGroup is \
            invalid. It must be None, a scalar, or a list/array or length \
            ncon=%d.' %(nCon))

        if upper is None:
            upper = [None for i in range(self.ncon)]
        elif numpy.isscalar(upper):
            upper = upper*numpy.ones(self.ncon)
        elif len(upper) == self.ncon:
            pass  # Some iterable object
        else:
            raise Error('The \'upper\' argument to addCon or addConGroup is \
            invalid. It must be None, a scalar, or a list/array or length \
            ncon=%d.' %(nCon))

        # ------ Process the scale argument
        scale = numpy.atleast_1d(scale)
        if len(scale) == 1:
            scale = scale[0]*numpy.ones(nCon)
        elif len(scale) == nCon:
            pass
        else:
            raise Error('The length of the \'scale\' argument to \
            addCon or addConGroup is %d, but the number of constraints \
            is %d.'%(len(scale), nCon))

        # Save lower and upper...they are only used for printing however
        self.lower = lower
        self.upper = upper
        
        # Now we determine what kind of constraint this is: 
        # 1. An equality constraint
        # 2. A upper bound on a 1-sided constraint
        # 3. A lower bound on a 1-sided constraint
        # 4. Lower and Upper bounds on 2-sided constraint
        # 5. No lower or upper bounds. Typically will only be used for
        # dummy constraint on an unconstrained problem.

        # The first 3, will give a single "constraint" in all
        # optimizers. Some optimizers can only do 1-sided constraints
        # so type 4 and 5 must be split into two separate constraints
        # automatically. 

        # This keeps track of the equality constraints:
        equalityConstraints = {'value': [], 'ind': [], 'fact': []}

        # All (inequality) constraints get added to
        # "twoSidedConstraints". This will be used in optimizers that
        # can do two-sided constraints properly
        twoSidedConstraints = {'lower': [], 'upper': [], 'ind': [], 'fact': []}

        # All (inequality) constraints are also added to
        # "oneSidedConstraints". These are processed such that the
        # lower bound is ALWAYS -INFINITY such that: con <= upper For
        # optimizers that need things <= zero, this can be processed
        # with a (-value) offset. One sided constraints need a fact
        # defined which is precisely 1.0 or -1.0. The -1.0 appears
        # when a greater-than-constraint is turned into a
        # less-than-constraint. 
        oneSidedConstraints = {'lower': [], 'upper': [], 'ind': [], 'fact': []}
        
        for icon in range(self.ncon):
            # Check for equality constraint:
            if lower[icon] == upper[icon] and lower[icon] is not None:
                equalityConstraints['value'].append(lower[icon]*scale[icon])
                equalityConstraints['ind'].append(icon)
                equalityConstraints['fact'].append(1.0)
                
            # Two sided constraint:
            elif lower[icon] is not None and upper[icon] is not None:
                twoSidedConstraints['lower'].append(lower[icon]*scale[icon])
                twoSidedConstraints['upper'].append(upper[icon]*scale[icon])
                twoSidedConstraints['ind'].append(icon)
                twoSidedConstraints['fact'].append(1.0)
                
                # TWO sets of 1 sided constraints:
                oneSidedConstraints['lower'].append(-INFINITY)
                oneSidedConstraints['upper'].append(upper[icon]*scale[icon])
                oneSidedConstraints['ind'].append(icon)
                oneSidedConstraints['fact'].append(1.0)

                oneSidedConstraints['lower'].append(-INFINITY)
                oneSidedConstraints['upper'].append(-lower[icon]*scale[icon])
                oneSidedConstraints['ind'].append(icon)
                oneSidedConstraints['fact'].append(-1.0)
                
            # Upper bound only:
            elif upper[icon] is not None:
                twoSidedConstraints['lower'].append(-INFINITY)
                twoSidedConstraints['upper'].append(upper[icon]*scale[icon])
                twoSidedConstraints['ind'].append(icon)
                twoSidedConstraints['fact'].append(1.0)
                
                # Just one, 1-sided constraint
                oneSidedConstraints['lower'].append(-INFINITY)
                oneSidedConstraints['upper'].append(upper[icon]*scale[icon])
                oneSidedConstraints['ind'].append(icon)
                oneSidedConstraints['fact'].append(1.0)

            # Lower bound only:
            elif lower[icon] is not None:
                twoSidedConstraints['lower'].append(lower[icon]*scale[icon])
                twoSidedConstraints['upper'].append(INFINITY)
                twoSidedConstraints['ind'].append(icon)
                twoSidedConstraints['fact'].append(1.0)
                
                # Just one, 1-sided constraint
                oneSidedConstraints['lower'].append(-INFINITY)
                oneSidedConstraints['upper'].append(-lower[icon]*scale[icon])
                oneSidedConstraints['ind'].append(icon)
                oneSidedConstraints['fact'].append(-1.0)

            # Fully unconstrained!
            elif lower[icon] is None and upper[icon] is None:
                twoSidedConstraints['lower'].append(-INFINITY)
                twoSidedConstraints['upper'].append(INFINITY)
                twoSidedConstraints['ind'].append(icon)
                twoSidedConstraints['fact'].append(1.0)
                                
                # Since this is just a dummy constraint, we only need
                # a single one....it can just be less than INFINITY
                oneSidedConstraints['lower'].append(-INFINITY)
                oneSidedConstraints['upper'].append(INFINITY)
                oneSidedConstraints['ind'].append(icon)
                oneSidedConstraints['fact'].append(1.0)
            # end if (con type)
        # end for (con loop)

        # Convert the stuff to arrays:
        oneSidedConstraints['ind'] = numpy.array(oneSidedConstraints['ind'],'intc')
        twoSidedConstraints['ind'] = numpy.array(twoSidedConstraints['ind'],'intc')
        equalityConstraints['ind'] = numpy.array(equalityConstraints['ind'],'intc')

        oneSidedConstraints['fact'] = numpy.array(oneSidedConstraints['fact'])
        twoSidedConstraints['fact'] = numpy.array(twoSidedConstraints['fact'])
        equalityConstraints['fact'] = numpy.array(equalityConstraints['fact'])

        equalityConstraints['value'] = numpy.array(equalityConstraints['value'])
        
        # Now save this information:
        self.equalityConstraints = equalityConstraints
        self.oneSidedConstraints = oneSidedConstraints
        self.twoSidedConstraints = twoSidedConstraints
        
    def __str__(self):
        """
        Print Constraint
        
        Documentation last updated:  April. 30, 2008 - Peter W. Jansen
        """
        res = ''
        for i in range(self.ncon):
            lower = self.lower[i]
            upper = self.upper[i]
            if lower is None:
                lower = 1e-20
            if upper is None:
                upper = 1e20
                
            res += '	 '+str(self.name).center(9) + \
                   '	  i %15.2e <= %8f <= %8.2e\n' %(
                lower, 0.0, upper)
       
        return res
