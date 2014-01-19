#!/usr/bin/env python
'''
pyOptSparse_optimization

Holds the Python Design Optimization Classes (base and inherited).

The class performs the basic functionality of
pyOpt/pyOpt_optimization, but eliminates most of the unused
functionality. The main purpose, therefore is to describe the
structure and sparsity pattery of a sparse optimization problem. 

Copyright (c) 2008-2013 by Dr. Gaetan Kenway
All rights reserved.

Developers:
-----------
- Dr. Gaetan K.W. Kenway (GKK)

History
-------
    v. 1.0  - Initial Class Creation (GKK, 2013)
'''

# =============================================================================
# Standard Python modules
# =============================================================================
import os, sys, time, copy
from collections import OrderedDict

# =============================================================================
# External Python modules
# =============================================================================
import numpy
from scipy import sparse

# =============================================================================
# Extension modules
# =============================================================================
from pyOpt_variable import Variable
from pyOpt_objective import Objective
from pyOpt_constraint import Constraint
from pyOpt_error import Error

# =============================================================================
# Misc Definitions
# =============================================================================
inf = 1e20  # define a value for infinity

# =============================================================================
# Optimization Class
# =============================================================================
class Optimization(object):
    '''
    Create a description of an optimization probelem. 

    Parameters
    ----------
    name : str
        Name given to optimization problem. This is name is currently
        not used for anything, but may be in the future.

    objFun : python function
        Python function handle of function used to evaluate the objective
        function.

    useGroups : bool
        Flag to specify whether or not design variables are returned in a
        flattened array or as a dictionary. It is **highly** recommened that
        useGroups is **always** used, even for small problems.
        '''
      
    def __init__(self, name, objFun, useGroups=True):

        self.name = name
        self.objFun = objFun
        self.useGroups = useGroups
        
        # Ordered dictionaries to keep track of variables and constraints
        self.variables = OrderedDict()
        self.constraints = OrderedDict()
        self.objectives = OrderedDict()
        self.dvOffset =  OrderedDict()

        # Sets to keep track of user-supplied names --- Keep track of
        # varSets, varGroup and constraint names independencely
        self.varSetNames = set()
        self.varGroupNames = set()
        self.conGroupNames = set()

        # Flag to determine if adding variables is legal. 
        self.ableToAddVariables = True

        return

    def _checkOkToAddVariables(self):
        if not self.ableToAddVariables:
            raise Error('No more variables can be added at this time. \
All variables must be added before constraints can be added.')
    
    def _finalizeDesignVariables(self):
        '''
        This function computes the ordering of the design variables
        and set the flag to prevent any more variables from being
        added '''

        dvCounter = 0

        for dvSet in self.variables.keys():
            # Check that varSet *actually* has variables in it:
            if len(self.variables[dvSet]) > 0:
                self.dvOffset[dvSet] = OrderedDict()
                self.dvOffset[dvSet]['n'] = [dvCounter, -1]
                for dvGroup in self.variables[dvSet]:
                    n = len(self.variables[dvSet][dvGroup])
                    self.dvOffset[dvSet][dvGroup] = [dvCounter, dvCounter + n, self.variables[dvSet][dvGroup][0].scalar]
                    dvCounter += n
                self.dvOffset[dvSet]['n'][1] = dvCounter
            else:
                # Get rid of the dvSet since it has no variable groups
                self.variables.pop(dvSet)
            # end if
        # end for
        self.ndvs = dvCounter
        self.ableToAddVariables = False

        return

    def addVarSet(self, name):
        '''An outer grouping of design variables. These sets are used
        when specifiying the sparsity structure of the constraint
        jacobian
        '''

        self._checkOkToAddVariables()

        if name in self.varSetNames:
            raise Error('The supplied name \'%s\' for a variable set \
            has already been used.'%(name))

        # end if
        self.varSetNames.add(name)
        self.variables[name]=OrderedDict()

        return

    def addVar(self, name, *args, **kwargs):
        '''
        This is a convience function. See addVarGroup for the
        appropriate list of parameters.
        '''

        self.addVarGroup(name, 1, *args, scalar=True, **kwargs)

        return 

    def addVarGroup(self, name, nVars, type='c', value=0.0, lower=None, upper=None, scale=1.0,
                    varSet='default', choices=None, **kwargs):
        '''
        Add a group of variables into a variable set. This is the main
        function used for adding variables to pyOptSparse.

        Parameters
        ----------
        name : str
            Name of variable group. This name should be unique across all the design variable groups

        nVars : int
            Number of design variables in this group.

        type : str. 
            String representing the type of variable. Suitable values for type
            are: 'c' for continuous variables, 'i' for integer values and
            'd' for discrete selection. 

        value : scalar or array. 
            Starting value for design variables. If it is a a scalar, the same
            value is applied to all 'nVars' variables. Otherwise, it must be
            iterable object with length equal to 'nVars'.

        lower : scalar or array. 
            Lower bound of variables. Scalar/array usage is the same as value
            keyword

        upper : scalar or array. 
            Upper bound of variables. Scalar/array usage is the same as value
            keyword
            
        scale : scalar or array. 
            Define a user supplied scaling variable for the design variable group.
            This is often necessary when design variables of widely varraying magnitudes
            are used within the same optimization. Scalar/array usage is the same
            as value keyword.

        varSet : str. 
            Specify which variable set this design variable group belons. If
            the variable set has not already been used, it will be added
            automatically. 

        choices : list
            Specify a list of choices for discrete design variables
            
        Examples
        --------
        >>> # Add a single design variable 'alpha' to the default variable set
        >>> optProb.addVar('alpha', type='c', value=2.0, lower=0.0, upper=10.0, \
        scale=0.1)
        >>> # Add a single variable to its own varSet
        >>> optProb.addVar('alpha_c1', type='c', value=2.0, lower=0.0, upper=10.0, \
        scale=0.1, varSet='alpha_c1')
        >>> # Add 10 unscaled variables of 0.5 between 0 and 1
        >>> optProb.addVarGroup('y', type='c', value=0.5, lower=0.0, upper=1.0, \
        scale=1.0, varSet='y_vars')
        >>> # Add another scaled variable to the varSet 'y_vars'
        >>> optProb.addVar('y2', type='c', value=0.25, lower=0.0, upper=1.0, \
        scale=.5, varSet='y_vars')
        
        Notes
        -----
        Calling addVar() and addVarGroup(..., nVars=1, ...) are
        **NOT** equilivant! The variable added with addVar() will be
        returned as scalar, while variable returned from addVarGroup
        will be an array of length 1.

        It is recommended that the addVar() and addVarGroup() calls
        follow the examples above by including all the keyword
        arguments. This make it very clear the itent of the script\'s
        author. The type, value, lower, upper and scale should be
        given for all variables even if the default value is used. 
        '''

        self._checkOkToAddVariables()

        if name in self.varGroupNames:
            raise Error('The supplied name \'%s\' for a variable group \
has already been used.'%(name))
        else:
            self.varGroupNames.add(name)
        # end if

        if not varSet in self.variables:
            self.addVarSet(varSet)
        # end if

        # Check that the type is ok:
        assert type in ['c','i','d'], 'Type must be one of \'c\' for continuous, \
\'i\' for integer or \'d\' for discrete.'
                    
        # ------ Process the value arguement
        value = numpy.atleast_1d(value)
        if len(value) == 1:
            value = value[0]*numpy.ones(nVars)
        elif len(value) == nVars:
            pass
        else:
            raise Error('The length of the \'value\' argument to \
 addVarGroup is %d, but the number of variables in nVars is %d.'%(len(value), nVars))
        # end if

        # ------ Process the lower bound argument
        if lower is None:
            lower = -inf*numpy.ones(nVars)
        else:
            lower = numpy.atleast_1d(lower)
            if len(lower) == 1:
                lower = lower[0]*numpy.ones(nVars)
            elif len(lower) == nVars:
                pass
            else:
                raise Error('The length of the \'lower\' argument to \
addVarGroup is %d, but the number of variables in nVars is %d.'%(len(lower), nVars))
            # end if
        # end if

        # ------ Process the upper bound argument
        if upper is None:
            upper = inf*numpy.ones(nVars)
        else:
            upper = numpy.atleast_1d(upper)
            if len(upper) == 1:
                upper = upper[0]*numpy.ones(nVars)
            elif len(upper) == nVars:
                pass
            else:
                raise Error('The length of the \'upper\' argument to \
addVarGroup is %d, but the number of variables in nVars is %d.'%(len(upper), nVars))
            # end if
        # end if

        # ------ Process the scale bound argument
        if scale is None:
            scale = numpy.ones(nVars)
        else:
            scale = numpy.atleast_1d(scale)
            if len(scale) == 1:
                scale = scale[0]*numpy.ones(nVars)
            elif len(scale) == nVars:
                pass
            else:
                raise Error('The length of the \'scale\' argument to \
addVarGroup is %d, but the number of variables in nVars is %d.'%(len(scale), nVars))
            # end if
        # end if

        # Determine if scalar i.e. it was called from addVar():
        scalar = kwargs.pop('scalar', False)

        # Now create all the variable objects
        self.variables[varSet][name] = []
        for iVar in xrange(nVars):
            varName = name + '_%d'%(iVar)
            self.variables[varSet][name].append(
                Variable(varName, type=type, value=value[iVar], lower=lower[iVar],
                         upper=upper[iVar], scale=scale[iVar], scalar=scalar, choices=choices))
        # end for

    def delVar(self, name):
        '''
        Delete a variable or variable group

        Parameters
        ----------
        name : str
           Name of variable or variable group to remove
           '''
        deleted = False
        for dvSet in self.variables.keys():
            for dvGroup in self.variables[dvSet]:
                if dvGroup == name:
                    self.variables[dvSet].pop(dvGroup)
                    deleted = True

        if not deleted:
            print '%s was not a valid design variable name'
            
    def delVarSet(self, name):
        '''
        Delete all variables belonging to a variable set

        Parameters
        ----------
        name : str
           Name of variable or variable group to remove
           '''

        assert name in self.variables, '%s not a valid varSet.'%name
        self.variables.pop(name)

    def addObj(self, name, *args, **kwargs):
        '''
        Add Objective into Objectives Set
        '''
        
        self.objectives[name] = Objective(name, *args, **kwargs)

        return
                
    def addCon(self, name, *args, **kwargs):
        '''
        Convenience function. See addConGroup() for more information
        '''
        
        self.addConGroup(name, 1, *args, **kwargs)

        return
                
    def addConGroup(self, name, nCon, lower=None, upper=None, scale=1.0,
                    linear=False, wrt=None, jac=None):
        '''
        Add a group of variables into a variable set. This is the main
        function used for adding variables to pyOptSparse.

        Parameters
        ----------

        name : str
            Constraint name. All names given to constraints must be unique

        nCon : int
            The number of constraints in this group

        lower : scalar or array
            The lower bound(s) for the constraint. If it is a scalar,
            it is applied to all nCon constraints. If it is an array,
            the array must be the same length as nCon.

        upper : scalar or array
            The upper bound(s) for the constraint. If it is a scalar,
            it is applied to all nCon constraints. If it is an array,
            the array must be the same length as nCon.

        scale : scalar or array

            A scaling factor for the constraint. It is generally
            advisible to have most optimization constraint around the
            same order of magnitude.

        linear : bool
            Flag to specifiy if this constraint is linear. If the
            constraint is linear, both the 'wrt' and 'jac' keyword
            arguments must be given to specify the constant portion of
            the constraint jacobian.

        wrt : iterable (list, set, OrderedDict, array etc)
            'wrt' stand for stands for 'With Respect To'. This
            specifies for what dvSets have non-zero jacobian values
            for this set of constraints. The order is not important.

        jac : dictionary
            For linear and sparse non-linear constraints, the constraint
            jacobian must be passed in. The structure is jac dictionary
            is as follows:

            {'dvSet1':<matrix1>, 'dvSet2', <matrix1>}

            They keys of the jacobian must correpsond to the dvSets
            givn in the wrt keyword argument. The dimensions of each
            "chunk" of the constraint jacobian must be consistent. For
            example, <matrix1> must have a shape of (nCon, nDvs) where
            nDVs is the **total** number of all design variables in
            dvSet1. <matrix1> may be a desnse numpy array or it may be
            scipy sparse matrix. It each case, the matrix shape must
            be as previously described. 

            Note that for nonlinear constraints (linear=False), the
            values themselves in the matrices in jac do not matter,
            but the sparsity structure **does** matter. It is
            imparative that entries that will at some point have
            non-zero entries have non-zero entries in jac
            argument. That is, we do not let the sparsity structure of
            the jacobian change throughout the optimization. This
            stipulation is automatically checked internally. 
            
              '''

        # If this is the first constraint, finalize the variables to
        # ensure no more variables can be added. 
        if self.ableToAddVariables:
            self._finalizeDesignVariables()

        if name in self.conGroupNames:
            raise Error('The supplied name \'%s\' for a constraint group \
has already been used.'%(name))

        self.conGroupNames.add(name)

        # ------ Process the lower bound argument
        if lower is None:
            lower = -inf*numpy.ones(nCon)
        else:
            lower = numpy.atleast_1d(lower)
            if len(lower) == 1:
                lower = lower[0]*numpy.ones(nCon)
            elif len(lower) == nCon:
                pass
            else:
                raise Error('The length of the \'lower\' argument to \
addConGroup is %d, but the number of constraints is %d.'%(len(lower), nCon))
            # end if
        # end if

        # ------ Process the upper bound argument
        if upper is None:
            upper = inf*numpy.ones(nCon)
        else:
            upper = numpy.atleast_1d(upper)
            if len(upper) == 1:
                upper = upper[0]*numpy.ones(nCon)
            elif len(upper) == nCon:
                pass
            else:
                raise Error('The length of the \'upper\' argument to \
addConGroup is %d, but the number of constraints is %d.'%(len(upper), nCon))
            # end if
        # end if

        # ------ Process the scale argument
        scale = numpy.atleast_1d(scale)
        if len(scale) == 1:
            scale = scale[0]*numpy.ones(nCon)
        elif len(scale) == nCon:
            pass
        else:
            raise Error('The length of the \'scale\' argument to \
 addConGroup is %d, but the number of constraints is %d.'%(len(scale), nCon))
        # end if
            
        # First check if 'wrt' is supplied...if not we just take all
        # the dvSet
        if wrt is None:
            wrt = self.variables.keys()
        else:
            # Sanitize the wrt input:
            if isinstance(wrt, str):
                wrt = [wrt.lower()]
            else: 
                try:
                    wrt = list(wrt)
                except:
                    raise Error('\'wrt\' must be a iterable list')
                # end try
            # end if
                    
            # Now, make sure that each dvSet the user supplied list
            # *actually* are DVsets
            for dvSet in wrt:
                if not dvSet in self.variables.keys():
                    raise Error('The supplied dvSet \'%s\' in \'wrt\' for the %s constraint, \
does not exist. It must be added with a call to addVar() or addVarGroup() with a dvSet=\'%s\' keyword argument.'%(dvSet, name, dvSet))
                # end if
            # end for
        # end if

        # Last thing for wrt is to reorder them such that dvsets are
        # in order. This way when the jacobian is assembled in
        # processDerivatives() the coorindate matrix will in the right
        # order.
        dvStart = []
        for dvSet in wrt:
            dvStart.append(self.dvOffset[dvSet]['n'][0])

        # This sort wrt using the keys in dvOffset
        wrt = [x for (y,x) in sorted(zip(dvStart, wrt))]

        # Now we know which DVsets this constraint will have a
        # derivative with respect to (i.e. what is in the wrt list)
            
        # Now, it is possible that jacobians were given for none, some
        # or all the dvSets defined in wrt. 
        if jac is None:

            # If the constraint is linear we have to *Force* the user to
            # supply a constraint jacobian for *each* of the values in
            # wrt. Otherwise, a matrix of zeros isn't meaningful for the
            # sparse constraints.

            if linear:
                raise Error('The \'jac\' keyword argument to addConGroup() must be supplied for a linear constraint')

            # without any additional information about the jacobian
            # structure, we must assume they are all dense. 
            jac = {}
            for dvSet in wrt:
                ss = self.dvOffset[dvSet]['n']                 
                ndvs = ss[1]-ss[0]
                jac[dvSet] = sparse.csr_matrix(numpy.ones((nCon, ndvs)))
                jac[dvSet].data[:] = 0.0
            # end for
                
            # Set a flag for the constraint object, that not returning them all is ok. 
            partialReturnOk = True

        else:
            # First sanitize input:
            if not isinstance(jac, dict):
                raise Error('The \'jac\' keyword argument to addConGroup() must be a dictionary')

            # Now loop over the set we *know* we need and see if any
            # are in jac. We will actually pop them out, and that way
            # if there is anything left at the end, we can tell the
            # user supplied information was unused. 
            tmp = copy.deepcopy(jac)
            jac = {}
            for dvSet in wrt:
                ss = self.dvOffset[dvSet]['n']                 
                ndvs = ss[1]-ss[0]

                try:
                    jac[dvSet] = tmp.pop(dvSet)
                    # Check that this user-supplied jacobian is in fact the right size
                except:
                    # No big deal, just make a dense component...and set to zero
                    jac[dvSet] = sparse.csr_matrix(numpy.ones((nCon, ndvs)))
                    jac[dvSet].data[:] = 0.0
                # end try
                    
                if jac[dvSet].shape[0] <> nCon or jac[dvSet].shape[1] <> ndvs:
                    raise Error('The supplied jacobian for dvSet \'%s\' in constraint %s, was the incorrect size. Expecting a jacobian\
 of size (%d,%d) but received a jacobian of size (%d,%d).'%(dvSet, name, nCon, ndvs, jac[dvSet].shape[0], jac[dvSet].shape[1]))
                # end if

                # Now check that the supplied jacobian is sparse of not:
                if sparse.issparse(jac[dvSet]):
                    # Excellent, the user supplied a sparse matrix or
                    # we just created one above. Convert to csr format
                    # if not already in that format.
                    jac[dvSet] = jac[dvSet].tocsr()
                else:
                    # Supplied jacobian is dense, replace any zero,
                    # before converting to csr format
                    jac[dvSet][numpy.where(jac[dvSet]==0)] = 1e-50
                    jac[dvSet] = sparse.csr_matrix(jac[dvSet])
                # end if
            # end for

            # If there is anything left in jac print a warning:
            for dvSet in tmp:
                print 'pyOptSparse Warning: An unused jacobian with dvSet key of \'%s\'\
 was unused. This will be ignored'%(dvSet)
            # end for

            # Finally partial returns NOT ok, since the user has
            # supplied information about the sparsity:
            partialReturnOk = False

        # end if

        # Scale the rows of each jacobian part:
        for dvSet in jac:
            self._csrRowScale(jac[dvSet], scale)

        # Finally! Create constraint object
        self.constraints[name] = Constraint(name, linear, wrt, jac, partialReturnOk,
                                            lower*scale, upper*scale, scale)

        return 
        
    def reorderConstraintJacobian(self, reorder=['nonLinear','linear']):
        '''
        ** This function should not need to be called by the end
           user**

        Reorder the supplied constraints such that all the nonlinear
        and linear constraints are grouped together. This function
        must be called by all optimizers regardless, since the con.rs
        and con.re values are computed in this function. 

        Parameters
        ----------
        reorder : list
            How to reorder:
            ['nonLinear','linear'] => put nonlinear ones first
            ['linear','nonLinear'] => put nonlinear ones first
            [] => Anything else: keep the same order as supplied.
      
        '''

        # Determine the total number of linear and nonlinear constraints:
        nlcon = 0 # Linear 
        nncon = 0 # nonlinear

        for iCon in self.constraints:
            if self.constraints[iCon].linear:
                nlcon += self.constraints[iCon].ncon
            else:
                nncon += self.constraints[iCon].ncon
            # end if
        # end for
        
        # Store number of linear and nonlinear constriants:
        self.nnCon = nncon
        self.nlCon = nlcon
        self.nCon = nncon + nlcon

        conScaleNonLinear = []
        conScaleFull = []

        # Loop over the constraints assigning the column start (cs)
        # and column end (ce) values. The actual ordering depends on
        # if constraints are reordered or not. 
        rowCounter = 0 
        if reorder == ['nonLinear','linear']:
            for iCon in self.constraints:
                con = self.constraints[iCon]
                if not con.linear:
                    con.rs = rowCounter
                    rowCounter += con.ncon
                    con.re = rowCounter
                    conScaleNonLinear.extend(con.scale)
                    conScaleFull.extend(con.scale)

            for iCon in self.constraints:
                con = self.constraints[iCon]
                if con.linear:
                    con.rs = rowCounter
                    rowCounter += con.ncon
                    con.re = rowCounter
                    conScaleFull.extend(con.scale)
        elif reorder == ['linear','nonLinear']:
            for iCon in self.constraints:
                con = self.constraints[iCon]
                if con.linear:
                    con.rs = rowCounter
                    rowCounter += con.ncon
                    con.re = rowCounter
                    conScaleFull.extend(con.scale)

            for iCon in self.constraints:
                con = self.constraints[iCon]
                if not con.linear:
                    con.rs = rowCounter
                    rowCounter += con.ncon
                    con.re = rowCounter
                    conScaleNonLinear.extend(con.scale)
                    conScaleFull.extend(con.scale)

        else: # No re-ordering
            for iCon in self.constraints:
                con = self.constraints[iCon]
                con.rs = rowCounter
                rowCounter += con.ncon
                con.re = rowCounter 
                if not con.linear:
                    conScaleNonLinear.extend(con.scale)
                conScaleFull.extend(con.scale)
            # end for
        # end if

        # Save constraint scaling arrays
        self.conScaleFull = numpy.array(conScaleFull)
        self.conScaleNonLinear = numpy.array(conScaleNonLinear)

        # Assemble the design variable scaling
        xscale = []
        for dvSet in self.variables.keys():
            for dvGroup in self.variables[dvSet]:
                for var in self.variables[dvSet][dvGroup]:
                    xscale.append(var.scale)
        self.xscale = numpy.array(xscale)

        return

    def processDerivatives(self, gobj_in, gcon_in, linearConstraints=False, 
                           nonlinearConstraints=True):
        '''
        ** This function should not need to be called by the end
        user**
        
        This generic function is used to assemble the objective
        gradient and the constraint jacobian.

        The two input flags are used to determine which if
        linear/nonlinear or both are included. Note that all cases the
        size of the jacobian is still (ncon x ndvs), ie the full
        size. However, only the requested entries (linear/nonlinear)
        are included. Also note that this function performs the pyOpt
        controlled scaling that is transparent to the user.

        Parameters
        ----------

        gobj_in : array or dict
            Objective gradient. Either a complete array or a dictionary of
            gradients given with respect to the dvSets

        gcon_in : array or dict
            Constraint gradients. Either a complete 2D array or a nested
            dictionary of gradients given with respect to the dvSets
        '''

        timeA = time.time()
        # Process the objective gradient. This may be vector or it may be
        # given as a dictionary. 
        dvSets = self.variables.keys()
        if isinstance(gobj_in, dict):
            gobj = numpy.zeros(self.ndvs)
            for key in gobj_in:
                # Check that the key matches something in dvSets
                if key in dvSets:
                    # Now check that the array is the correct length:

                    # ss = start/stop is a length 2 array of the
                    # indices for dvSet given by key
                    ss = self.dvOffset[key]['n'] 
                    if len(gobj_in[key]) == ss[1]-ss[0]:
                        # Everything checks out so set:
                        gobj[ss[0]:ss[1]] = gobj_in[key]
                    else:
                        raise Error('The length of the objective deritative for dvSet %s\
 is the incorrect length. Expecting a length of %d but received a length of %d.'%(dvSet, ss[1]-ss[0], len(gobj_in[key])))
                    # end if
                else:
                    print 'Warning: The key \'%s\' in g_obj does not match any of the added DVsets. \
This derivative will be ignored'%(key)
                # end if
            # end for
        else:
            # Otherwise we will assume it is a vector:
            gobj = numpy.atleast_1d(gobj_in).copy()
            if len(gobj) <> self.ndvs:
                raise Error('The length of the objective derivative for all design variables\
 is not the correct size. Received size %d, should be size %d.'%(len(gobj), self.ndvs))
        # end if
                
        # Finally scale the objective gradient based on the scaling data.
        gobj *= self.xscale

        # If the user has supplied a complete dense numpy array for
        # the jacobain AND all the constriants are dense
        if isinstance(gcon_in, numpy.ndarray):
            if gcon_in.shape == (self.nCon, self.ndvs):
                gcon_in[numpy.where(gcon_in==0)] = 1e-50
                # Don't forget to scale:
                tmp = gcon_in.copy()
                for i in xrange(self.ndvs):
                    tmp[:,i] *= self.xscale[i]

                # Do constraint scaling and convert to coo
                gcon = sparse.coo_matrix(tmp)
                self._cooRowScale(gcon, self.conScaleFull)
                print 'Processing Derviatives Time: %8.3f seconds.'%(time.time()-timeA)

                # Quick Return here since we're done
                return gobj, gcon

            else:
                raise Error('The dense jacobian return was the incorrect size. Expecting \
size of (%d, %d) but received size of (%d, %d).'%(self.nCon, self.ndvs, gcon_in.shape[0], gcon_in.shape[1]))
            # end if
        # end if

        # Data for storing the values in COOrdinate format
        data = []
        row  = []
        col  = []

        # Otherwise, process constraints in the dictionary form. 
        # Loop over all constraints:
        for iCon in self.constraints:
            con = self.constraints[iCon]
            if ((linearConstraints and con.linear) or
                (nonlinearConstraints and not con.linear)):

                if not con.name in gcon_in:
                    raise Error('The jacobian for the constraint \'%s\' was \
not found in the returned dictionary.'%con.name)
                # end if

                if not con.partialReturnOk:
                    # The keys in gcon_in[iCon] MUST match PRECISELY
                    # the keys in con.wrt....The user told us they
                    # would supply derivatives wrt to these sets, and
                    # then didn't, so scold them. 
                    for dvSet in con.jac.keys():
                        if dvSet not in gcon_in[iCon]:
                            raise Error('Constraint \'%s\' was expecting a jacobain with respect to dvSet \'%s\' as \
was supplied in addConGroup(). This was not found in the constraint jacobian dictionary'%(con.name, dvSet))
                    # end for
                # end if

                # We assume gcon_in MUST be a dictionary
                for key in con.wrt:

                    ss = self.dvOffset[key]['n'] 
                    ndvs = ss[1]-ss[0]
                    if key in gcon_in[iCon]:
                        if sparse.issparse(gcon_in[iCon][key]):
                            # Excellent, the user supplied a sparse matrix
                            # Convert to csr format if not already in that
                            # format.
                            tmp = gcon_in[iCon][key].copy().tocsr()
                        else:
                            # Supplied jacobian is dense, replace any zero,
                            # before converting to csr format
                            gcon_in[iCon][key][numpy.where(gcon_in[iCon][key]==0)] = 1e-50
                            tmp = sparse.csr_matrix(gcon_in[iCon][key].copy())
                        # end if
                    else:
                        # Just use stored jacobian that contains just zeros:
                        tmp = con.jac[key]
                    # end if

                    # Now check that the jacobian is the correct shape
                    if not(tmp.shape[0] == con.ncon and tmp.shape[1] == ndvs):
                        raise Error('The shape of the supplied constraint jacobian for constraint %s is incorrect. \
 Expected an array of shape (%d,%d), but received an array of shape (%d, %d).'%(con.name, con.ncon, ndvs, tmp.shape[0], tmp.shape[1]))

                    # Now check that the csr matrix has the correct number of non zeros:
                    if tmp.nnz <> con.jac[key].nnz:
                        raise Error('The number of nonzero elements for \
  constraint group \'%s\' was not the correct size. The supplied jacobian has \
 %d nonzero entries, but must contain %d nonzero entries.'%(con.name, tmp.nnz, con.jac[key].nnz))
                    # end if

                    # Loop over rows in constraint jacobian:
                    for iRow in xrange(con.ncon):
                        # Loop over the number of nonzero entries in this row:
                        for ii in xrange(con.jac[key].indptr[iRow], con.jac[key].indptr[iRow+1]):
                            row.append(con.rs + iRow)
                            icol = self.dvOffset[key]['n'][0] + con.jac[key].indices[ii]
                            col.append(icol)
                            data.append(tmp.data[ii]*self.xscale[icol])
                        # end for
                    # end for
                # end for
            # end if
        # end for

        # Create coo matrix and scale the rows
        gcon = sparse.coo_matrix((data, (row, col)),(self.nCon, self.ndvs))
        self._cooRowScale(gcon, self.conScaleFull)

        print 'Processing Derviatives Time: %8.3f seconds.'%(time.time()-timeA)

        return gobj, gcon

    def processConstraints(self, fcon_in, linearConstraints=False, 
                           nonlinearConstraints=True):
        '''
        ** This function should not need to be called by the end
        user**

        Parameters
        ----------
        fcon_in : array or dict
            Array of constraint values or a dictionary of constraint
            values

        linearConstraints : bool
            Flag as to whether or not linear constraints are to be included

        nonLinearConstraint : bool
            Flag as to whether or not nonlinear constraints are to be included. 
            '''

        # We will actually be a little leniant here; the user CAN
        # return an iterable of the correct length and we will accept
        # that. Otherwise we will use the dictionary formulation

        if not isinstance(fcon_in, dict):
            fcon = numpy.atleast_1d(fcon_in)
            if len(fcon) == self.nnCon:
                return self.conScaleNonLinear*fcon
            else:
                raise Error('The constraint array was the incorrect size. \
It must contain %d elements (nonlinear constraints only), but an arrary of \
size %d was given.'%(self.nnCon, len(fcon)))
            # end if
        else:
            # Process as a dictionary:
            # Loop over (nonlinear) constraints and extract as required:
            fcon = []
            for iCon in self.constraints:
                if not self.constraints[iCon].linear:
                    if iCon in fcon_in:
                        # Make sure it is at least 1dimension:
                        c = numpy.atleast_1d(fcon_in[iCon])
                        
                        # Make sure it is the correct size:
                        if len(c) == self.constraints[iCon].ncon:
                            fcon.extend(c)
                        else:
                            raise Error('%d constraint values were returned \
    %s, but expected %d.'%(len(fcon_in[iCon]), iCon, self.variables[iCon].ncon))
                        # end if
                    else:
                        raise Error('No constraint values were found for the \
constraint %s.'%(iCon))
                    # end if
                # end if
            # end for

            # Finally convert to array and scale
            fcon = self.conScaleNonLinear*numpy.array(fcon)
        # end if
 
        return fcon

    def convertToDense(self):
        '''
        ** This function should not need to be called by the end
        user**
        
        Take a sparse optimization problem definition and convert to a
        dense representation. This function will be called
        automatically by optimizers that do not support sparse
        jacobians or do not support linear constraints.
        '''

        return

    def processX(self, x):
        '''
        ** This function should not need to be called by the end
        user**

        Take the flattened array of variables in \'x\' and return a
        dictionary of variables keyed on the name of each variable
        group if useGroups is True. 

        Parameters
        ----------
        x : array
            Flattened array from optimizer
        '''
        if self.useGroups:
            xg = {}
            for dvSet in self.variables.keys():
                for dvGroup in self.variables[dvSet]:
                    istart = self.dvOffset[dvSet][dvGroup][0]
                    iend   = self.dvOffset[dvSet][dvGroup][1]
                    scalar = self.dvOffset[dvSet][dvGroup][2]
                    xg[dvGroup] = x[istart:iend]
                    if scalar:
                        xg[dvGroup] = x[istart]
                # end for
            # end for
            return xg 
        else:
            return x
        # end if

    def _csrRowScale(self, mat, vec):
        '''
        Scale rows in csr matrix. Amazingly enough this is effectively
        impossible with scipy.sparse if you want to keep the nonzero
        structure. So we will brute force it here.
        '''
        assert mat.shape[0] == len(vec)
        for iRow in xrange(mat.shape[0]):
            mat.data[mat.indptr[iRow]:mat.indptr[iRow+1]] *= vec[iRow]
        # end for

    def _cooRowScale(self, mat, vec):
        ''' 
        Scale rows of coo matrx. See _csrRowScale for why
        '''
        assert mat.shape[0] == len(vec)
        for i in xrange(len(mat.data)):
            mat.data[i] *= vec[mat.row[i]]

    def __str__(self):
        '''
        Print Structured Optimization Problem
        '''
        
        text = '''\nOptimization Problem -- %s\n%s\n
        Objective Function: %s\n\n    Objectives:
        Name        Value        Optimum\n''' %(self.name,'='*80,self.objFun.__name__)

        for obj in self.objectives:
            lines = str(self.objectives[obj]).split('\n')
            text += lines[1] + '\n'
        # end for

        text += '''\n	Variables (c - continuous, i - integer, d - discrete):
        Name    Type       Value       Lower Bound  Upper Bound\n'''

        for dvSet in self.variables:
            for dvGroup in self.variables[dvSet]:
                for var in self.variables[dvSet][dvGroup]:
                    lines = str(var).split('\n')
                    text+= lines[1] + '\n'
                # end for
            # end for
        # end for

        print '	    Name        Type'+' '*25+'Bound\n'+'	 '
        if len(self.constraints) > 0:
            text += '''\n	Constraints (i - inequality, e - equality):
        Name    Type                    Bounds\n'''
            for iCon in self.constraints:
                text += str(self.constraints[iCon])
            # end for
        # end if
        
        return text
 
    def printSparsity(self):
        '''
        This function prints an (ascii) visualization of the jacobian
        sparsity structure. This helps the user visualize what
        pyOptSparse has been given and helps ensure it is what the
        user expected. It is highly recommended this function be
        called before the start of every optimization to verify the
        optimization problem setup.
        '''

        # Header describing what we are printing:
        print '+'+'-'*78+'-'+'+'
        print '|' + ' '*19 +'Sparsity structure of constraint Jacobian' + ' '*19 + '|'
        print '+'+'-'*78+'-'+'+'

        # We will do this with a 2d numpy array of characters since it
        # will make slicing easier

        # First determine the requried number of rows 
        nRow = 1 # Header
        nRow += 1 # Line
        maxConNameLen = 0
        hasLinear = False
        for iCon in self.constraints:
            nRow += 1 # Name
            con = self.constraints[iCon]
            maxConNameLen = max(maxConNameLen, len(con.name)+3+int(numpy.log10(con.ncon))+1)
            nRow += 1 # Line
            if self.constraints[iCon].linear:
                hasLinear = True
        # end for
        if hasLinear:
            nRow += 1 # Extra line to separate linear constraints

        # And now the columns:
        nCol = maxConNameLen
        nCol += 2 # Space plus line
        varCenters = []
        for iVar in self.variables:
            nvar = self.dvOffset[iVar]['n'] [1] - self.dvOffset[iVar]['n'][0]
            var_str = iVar + ' (%d)'%nvar

            varCenters.append(nCol + len(var_str)/2 + 1)
            nCol += len(var_str)
            nCol += 2 # Spaces on either side
            nCol += 1 # Line 

        # end for
        txt = numpy.zeros((nRow, nCol),dtype=str)
        txt[:,:] = ' '
        # Outline of the matrix on left and top
        txt[1,maxConNameLen+1:-1] = '-'
        txt[2:-1,maxConNameLen+1] = '|'
     
        # Print the variable names:
        iCol = maxConNameLen + 2
        for iVar in self.variables:
            nvar = self.dvOffset[iVar]['n'] [1] - self.dvOffset[iVar]['n'][0]
            var_str = iVar + ' (%d)'%nvar
            l = len(var_str)
            txt[0, iCol+1 :iCol + l+1] = list(var_str)
            txt[2:-1, iCol + l + 2] = '|'
            iCol += l + 3

        # Print the constraint names;
        iRow = 2

        # Do the nonlinear ones first:
        for iCon in self.constraints:
            con = self.constraints[iCon]
            if not con.linear:
                name = con.name + ' (%d)'%con.ncon
                l = len(name)
                # The name
                txt[iRow, maxConNameLen-l:maxConNameLen] = list(name)

                # Now we write a 'D' for dense, 'S' for sparse or nothing. 
                varKeys = self.variables.keys()
                for iVar in xrange(len(varKeys)):
                    if varKeys[iVar] in con.wrt:
                        txt[iRow, varCenters[iVar]] = 'X'

                # The separator
                txt[iRow+1, maxConNameLen+1:] = '-'
                iRow += 2
            # end if
        # end for
            
        # Print an extra '---' to distinguish:
        if hasLinear:
            txt[iRow, maxConNameLen+1:] = '-'
            iRow += 1

        # Do the nonlinear ones first and then the linear ones:
        for iCon in self.constraints:
            con = self.constraints[iCon]
            if con.linear:
                name = con.name + ' (%d)'%con.ncon
                l = len(name)
                # The name
                txt[iRow, maxConNameLen-l:maxConNameLen] = list(name)

                # Now we write a 'D' for dense, 'S' for sparse or nothing. 
                varKeys = self.variables.keys()
                for iVar in xrange(len(varKeys)):
                    if varKeys[iVar] in con.wrt:
                        txt[iRow, varCenters[iVar]] = 'X'
                    # end if
                # end for


                # The separator
                txt[iRow+1, maxConNameLen+1:] = '-'
                iRow += 2
            # end if
        # end for

        # Corners
        txt[1, maxConNameLen+1] = '+'
        txt[-1,maxConNameLen+1] = '+'
        txt[1,-1] = '+'
        txt[-1,-1] = '+'
        for i in xrange(len(txt)):
            print ''.join(txt[i])
        # end for

        return 

#==============================================================================
# Optimization Test
#==============================================================================
if __name__ == '__main__':
    
    print 'Testing Optimization...'
    optprob = Optimization('Optimization Problem',{})
    
    
