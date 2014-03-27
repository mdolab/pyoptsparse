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
import copy
import numpy
from scipy import sparse
from .pyOpt_error import Error
INFINITY = 1e20
eps = numpy.finfo(1.0).eps
# =============================================================================
# Constraint Class
# =============================================================================
class Constraint(object):
    """
    Constraint Class Initialization
    """
    def __init__(self, name, nCon, linear, wrt, jac, lower, upper, scale):

        self.name = name
        self.ncon = nCon
        self.linear = linear
        self.wrt = wrt
        self.jac = jac
        self.partialReturnOk = None
        self.scale = scale
        self.rs = None
        self.re = None
        # Before we can do the processing below we need to have lower
        # and upper arguments expanded:

        if lower is None:
            lower = [None for i in range(self.ncon)]
        elif numpy.isscalar(lower):
            lower = lower*numpy.ones(self.ncon)
        elif len(lower) == self.ncon:
            pass  # Some iterable object
        else:
            raise Error("The 'lower' argument to addCon or addConGroup is "
                        "invalid. It must be None, a scalar, or a "
                        "list/array or length ncon=%d." % nCon)

        if upper is None:
            upper = [None for i in range(self.ncon)]
        elif numpy.isscalar(upper):
            upper = upper*numpy.ones(self.ncon)
        elif len(upper) == self.ncon:
            pass  # Some iterable object
        else:
            raise Error("The 'upper' argument to addCon or addConGroup is "
                        "invalid. It must be None, a scalar, or a "
                        "list/array or length ncon=%d." %nCon)

        # ------ Process the scale argument
        scale = numpy.atleast_1d(scale)
        if len(scale) == 1:
            scale = scale[0]*numpy.ones(nCon)
        elif len(scale) == nCon:
            pass
        else:
            raise Error("The length of the 'scal\' argument to "
                        "addCon or addConGroup is %d, but the number of "
                        "constraints is %d."% (len(scale), nCon))

        # Save lower and upper...they are only used for printing however
        self.lower = lower
        self.upper = upper
        # The current value of the constraint (for printing purposes)
        self.value = numpy.zeros(self.ncon)
        
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

    def finialize(self, variables, dvOffset, index):
        """ **This function should not need to be called by the user**
        
        After the design variables have been finialized and the order
        is known we can check the constraint for consistency. 

        Parameters
        ----------
        variables : Ordered Dict
            The pyOpt variable list after they have been finialized.

        dvOffset : dict
            Design variable offsets from pyOpt_optimization

        index : int
            The starting index of this constraint in natural order

        """
        
        # Set the row start and end
        self.rs = index
        self.re = index + self.ncon

        # First check if 'wrt' is supplied...if not we just take all
        # the dvSet
        if self.wrt is None:
            self.wrt = list(variables.keys())
        else:
            # Sanitize the wrt input:
            if isinstance(self.wrt, str):
                self.wrt = [self.wrt.lower()]
            else: 
                try:
                    self.wrt = list(self.wrt)
                except:
                    raise Error("The 'wrt' argument to constraint '%s' must "
                                "be an iterable list"% self.name)

            # We allow 'None' to be in the list...they are null so
            # just pop them out:
            self.wrt = [dvSet for dvSet in self.wrt if dvSet != None]
                    
            # Now, make sure that each dvSet the user supplied list
            # *actually* are DVsets
            for dvSet in self.wrt:
                if not dvSet in variables:
                    raise Error("The supplied dvSet '%s' in 'wrt' "
                                "for the %s constraint, does not exist. It "
                                "must be added with a call to addVar() or "
                                "addVarGroup() with a dvSet='%s' keyword "
                                "argument."% (dvSet, self.name, dvSet))

        # Last thing for wrt is to reorder them such that dvsets are
        # in order. This way when the jacobian is assembled in
        # processDerivatives() the coorindate matrix will in the right
        # order.
        dvStart = []
        for dvSet in self.wrt:
            dvStart.append(dvOffset[dvSet]['n'][0])

        # This sort wrt using the keys in dvOffset
        self.wrt = [x for (y, x) in sorted(zip(dvStart, self.wrt))]

        # Now we know which DVsets this constraint will have a
        # derivative with respect to (i.e. what is in the wrt list)
            
        # Now, it is possible that jacobians were given for none, some
        # or all the dvSets defined in wrt. 
        if self.jac is None:
            # If the constraint is linear we have to *Force* the user to
            # supply a constraint jacobian for *each* of the values in
            # wrt. Otherwise, a matrix of zeros isn't meaningful for the
            # sparse constraints.

            if self.linear:
                raise Error("The 'jac' keyword to argument to addConGroup() "
                            "must be supplied for a linear constraint. "
                            "The constraint in error is %s."% self.name)

            # without any additional information about the jacobian
            # structure, we must assume they are all dense. 
            self.jac = {}
            for dvSet in self.wrt:
                ss = dvOffset[dvSet]['n']                 
                ndvs = ss[1]-ss[0]
                self.jac[dvSet] = sparse.coo_matrix(numpy.ones((self.ncon, ndvs)))
                self.jac[dvSet].data[:] = 1e-50
                
            # Set a flag for the constraint object, that not returning
            # them all is ok.
            self.partialReturnOk = True

        else:
            # First sanitize input:
            if not isinstance(self.jac, dict):
                raise Error("The 'jac' keyword argument to addConGroup() "
                            "must be a dictionary. The constraint in error "
                            "is %s."% self.name)

            # Now loop over the set we *know* we need and see if any
            # are in jac. We will actually pop them out, and that way
            # if there is anything left at the end, we can tell the
            # user supplied information was unused. 
            tmp = copy.deepcopy(self.jac)
            self.jac = {}
            for dvSet in self.wrt:
                ss = dvOffset[dvSet]['n']                 
                ndvs = ss[1]-ss[0]

                try:
                    self.jac[dvSet] = tmp.pop(dvSet)
                    # Check that this user-supplied jacobian is in
                    # fact the right size
                except:
                    # No big deal, just make a dense component...and
                    # set to zero
                    self.jac[dvSet] = sparse.coo_matrix(numpy.ones((self.ncon, ndvs)))
                    self.jac[dvSet].data[:] = 1e-50
                    
                if (self.jac[dvSet].shape[0] != self.ncon or 
                    self.jac[dvSet].shape[1] != ndvs):
                    raise Error("The supplied jacobian for dvSet \%s' "
                                "in constraint %s, was the incorrect size. "
                                "Expecting a jacobian of size (%d, %d) but "
                                "received a jacobian of size (%d, %d)."%(
                                    dvSet, self.name, self.ncon, ndvs, 
                                    self.jac[dvSet].shape[0],
                                    self.jac[dvSet].shape[1]))

                # Now check that the supplied jacobian is sparse of not:
                if sparse.issparse(self.jac[dvSet]):
                    # Excellent, the user supplied a sparse matrix or
                    # we just created one above. Convert to csr format
                    # if not already in that format.
                    self.jac[dvSet] = self.jac[dvSet].tocoo()
                else:
                    # Supplied jacobian is dense, replace any zero, 
                    # before converting to csr format
                    self.jac[dvSet][numpy.where(self.jac[dvSet]==0)] = 1e-50
                    self.jac[dvSet] = sparse.coo_matrix(self.jac[dvSet])
            # end for (dvSet)

            # If there is anything left in jac print a warning:
            for dvSet in tmp:
                print("pyOptSparse Warning: A jacobian with dvSet key of "
                      "'%s' was unused in constraint %s. This will be "
                      "ignored."% ( dvSet, self.name))

            # Since this function *may* be called multiple times, only
            # set paritalReturnOk if it was the first pass:
            if self.partialReturnOk is None:
                # Finally partial returns NOT ok, since the user has
                # supplied information about the sparsity:
                self.partialReturnOk = False

        # end if (if Jac)

    def __str__(self):
        """
        Print Constraint
        
        Documentation last updated:  April. 30, 2008 - Peter W. Jansen
        """
        res = ''
        for i in range(self.ncon):
            lower = self.lower[i]
            upper = self.upper[i]
            value = self.value[i]
            if lower is None:
                lower = 1e-20
            if upper is None:
                upper = 1e20
                
            res += '	 '+str(self.name).center(9) + \
                   '	  i %15.2e <= %8f <= %8.2e\n' %(
                lower, value, upper)
       
        return res
