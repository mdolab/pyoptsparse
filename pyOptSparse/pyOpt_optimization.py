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
Revision: 1.0   $Date: 19/09/2013 21:00$


Developers:
-----------
- Dr. Gaetan K.W. Kenway (GKK)

History
-------
    v. 1.0  - Initial Class Creation (GKK, 2013)

'''

__version__ = '$Revision: $'

'''
To Do:


'''

# =============================================================================
# Standard Python modules
# =============================================================================
import os, sys, time
try:
    from collections import OrderedDict
except ImportError:
    # python 2.6 or earlier, use backport
    # get using "pip install ordereddict"
    from ordereddict import OrderedDict
    

# =============================================================================
# External Python modules
# =============================================================================
import numpy
from scipy import sparse
# =============================================================================
# Extension modules
# =============================================================================
from pyOptSparse import Variable
from pyOptSparse import Objective
from pyOptSparse import Constraint

# =============================================================================
# Misc Definitions
# =============================================================================
inf = 10.E+20  # define a value for infinity

# =============================================================================
# Optimization Class
# =============================================================================
class Optimization(object):
    
    '''
    Optimization Problem Class
    '''
    
    def __init__(self, name, obj_fun, use_groups=False, *args, **kwargs):
        
        '''
        Optimization Problem Class Initialization
        
        **Arguments:**
        
        - name -> STR: Solution name
        - obj_fun -> FUNC: Objective function
        
        **Keyword arguments:**
        
        - var_set -> INST: Variable set, *Default* = None
        - obj_set -> INST: Objective set, *Default* = None
        - con_set -> INST: Constraints set, *Default* = None
        - use_groups -> BOOL: Use of group identifiers flag, *Default* = False
        
        Documentation last updated:  May. 23, 2011 - Ruben E. Perez
        '''
        
        # 
        self.name = name
        self.obj_fun = obj_fun
        self.use_groups = use_groups
        
        # Ordered dictionaries to keep track of variables and constraints
        self.variables = OrderedDict()
        self.constraints = OrderedDict()
        self.objectives = OrderedDict()
        self.dvOffset =  OrderedDict()

        # Separate (ordered) dictionaries to keep track of each of the
        # separate types of constraints: Dense (assumed nonlinear),
        # Sparse nonlinear, and sparseLinear
        self.denseConstraints = OrderedDict()
        self.sparseConstraints = OrderedDict()
        self.linearSparseConstraints = OrderedDict()

        # A set to keep track of user-supplied names --- Keep track of
        # varSets, varGroup and constraint names independencely
        self.varSetNames = set()
        self.varGroupNames = set()
        self.conGroupNames = set()

        # Flag to determine if adding variables is legal. 
        self.ableToAddVariables = True

        return

    def _checkOkToAddVariables(self):
        if not self.ableToAddVariables:
            print 'Error: No more variables can be added at this time. \
All variables must be added before constraints can be added.'
            sys.exit(1)
        # end if
        
        return
    
    def _finalizeDesignVariables(self):
        '''
        This function computes the ordering of the design variables
        and set the flag to prevent any more variables from being
        added '''

        dvCounter = 0
        for dvSet in self.variables.keys():
            self.dvOffset[dvSet] = OrderedDict()
            self.dvOffset[dvSet]['n'] = [dvCounter, -1]
            for dvGroup in self.variables[dvSet]:
                n = len(self.variables[dvSet][dvGroup])
                self.dvOffset[dvSet][dvGroup] = [dvCounter, dvCounter + n, self.variables[dvSet][dvGroup][0].scalar]
                dvCounter += n
            self.dvOffset[dvSet]['n'][1] = dvCounter
        # end for
        self.ndvs = dvCounter
        self.ableToAddVariables = False

        return

    def addVarSet(self, name):
        '''An outer grouping of design variables. These sets are used
        when specifiying the sparsity structuer of the constraint
        jacobian
        '''

        self._checkOkToAddVariables()
        
        if name in self.varSetNames:
            print 'Error: The supplied name \'%s\' for a variable set \
has already been used.'%(name)
            return
        # end if
        self.varSetNames.add(name)
        self.variables[name]=OrderedDict()

        return

    def addVar(self, name, *args, **kwargs):
        '''
        convenience function. See addVarGroup for more information
        '''

        self.addVarGroup(name, 1, *args, scalar=True, **kwargs)

        return 

    def addVarGroup(self, name, nvars, type='c', value=0.0, varSet='default', scalar=False, **kwargs):
        
        '''
        Add a Group of Variables into Variables Set
        
        **Arguments:**
        
        - name -> STR: Variable Group Name
        - nvars -> INT: Number of variables in group
        
        **Keyword arguments:**
        
        - type -> STR: Variable type ('c'-continuous, 'i'-integer, 'd'-discrete), *Default* = 'c'
        - value ->INT/FLOAT: Variable starting value, *Default* = 0.0
        
        Documentation last updated:  Feb. 07, 2011 - Peter W. Jansen
        '''

        self._checkOkToAddVariables()

        if name in self.varGroupNames:
            print 'Error: The supplied name \'%s\' for a variable group \
has already been used.'%(name)
            return
        else:
            self.varGroupNames.add(name)
        # end if

        if not varSet in self.variables:
            self.addVarSet(varSet)
        # end if
            
        # ------ Process the value arguement
        value = numpy.atleast_1d(value)
        if len(value) == 1:
            value = value[0]*numpy.ones(nvars)
        elif len(value) == nvars:
            pass
        else:
            print 'Error: The length of the \'value\' argument to \
 addVarGroup is %d, but the number of variables in nvars is %d.'%(len(value), nvars)
            sys.exit(1)
        # end if

        # ------ Process the lower bound argument
        lower = numpy.atleast_1d(kwargs.pop('lower', -inf*numpy.ones(nvars)))
        if len(lower) == 1:
            lower = lower[0]*numpy.ones(nvars)
        elif len(lower) == nvars:
            pass
        else:
            print 'Error: The length of the \'lower\' argument to \
 addVarGroup is %d, but the number of variables in nvars is %d.'%(len(lower), nvars)
            sys.exit(1)
        # end if

        # ------ Process the upper bound argument
        upper = numpy.atleast_1d(kwargs.pop('upper', inf*numpy.ones(nvars)))
        if len(upper) == 1:
            upper = upper[0]*numpy.ones(nvars)
        elif len(upper) == nvars:
            pass
        else:
            print 'Error: The length of the \'upper\' argument to \
 addVarGroup is %d, but the number of variables in nvars is %d.'%(len(upper), nvars)
            sys.exit(1)
        # end if

        # ------ Process the scale argument
        scale = numpy.atleast_1d(kwargs.pop('scale', numpy.ones(nvars)))
        if len(scale) == 1:
            scale = scale[0]*numpy.ones(nvars)
        elif len(scale) == nvars:
            pass
        else:
            print 'Error: The length of the \'scale\' argument to \
 addVarGroup is %d, but the number of variables in nvars is %d.'%(len(scale), nvars)
            sys.exit(1)
        # end if

        # Now create all the variable objects
        self.variables[varSet][name] = []
        for iVar in xrange(nvars):
            varName = name + '_%d'%(iVar)
            self.variables[varSet][name].append(
                Variable(varName, type=type, value=value[iVar]*scale[iVar], 
                         lower=lower[iVar]*scale[iVar], upper=upper[iVar]*scale[iVar]))
            self.variables[varSet][name][-1].scale = scale[iVar]
            self.variables[varSet][name][-1].scalar = scalar
        # end for
                                                
    def addObj(self, name, *args, **kwargs):
        
        '''
        Add Objective into Objectives Set
        
        Documentation last updated:  March. 27, 2008 - Ruben E. Perez
        '''
        
        self.objectives[name] = Objective(name, *args, **kwargs)

        return
                
    def addCon(self, name, *args, **kwargs):
        '''
        convenience function. See addConGroup for more information
        '''
        
        self.addConGroup(name, 1, *args, **kwargs)

        return
                
    def addConGroup(self, name, ncons, type='i', dense=True,
                    linear=False, wrt=None, jac=None, **kwargs):
        
        '''
        Add a Group of Constraints into Constraints Set
        
        **Arguments:**
        
        - name -> STR: Constraint group name
        - ncons -> INT: Number of constraints in group
        
        **Keyword arguments:**
        
        - type -> STR: Constraint type ('i'-inequality, 'e'-equality), *Default* = 'i'
        - Only inequality constriants are supported 
        
        - dense-> Boolean: Whether or not this constraint is taken to
          be dense or not. If the constraint it sparse the
          'wrt' and 'jac' keyword arguments must be provided

        - linear -> Boolean: Whether or not this constraint is
          linear. If it is linear and not dense, the jacobian provided
          must contain the actual (fixed) jacobian values

        - wrt -> iterable (list, set, OrderedDict, array etc): wrt
          stands for 'With Respect To'. This specifies for what dvSets
          have non zero jacobian values for this set of
          constraints. The order is important. For example, lets
          assume the user has defined dvSet 'A' and dvSet 'B'. 'A' has
          10 variables and 'B' has 15 variables and ncons (in this
          function) is 25. If the user calls:

          optProb.addConGroup('aSetOfConstraints', 25, type='i', dense=False,
          linear=False, wrt=['A','B'], jac=jac, lower=lower, upper=upper)

          the jacobian (jac) must have shape of (25, 15). It is
          expected the first 10 columns correspond to the derivates
          wrt to variables 'A' and the last 5 columns wrt variables
          'B'. If wrt was given as wrt=['B','A'] the opposite order
          would be used. Note that the 'wrt' order MUST be used in all
          subequent returns of the jacobian.

          jac -> scipy.sparse matrix: The sparse matrix jacobian. For
          nonlinear constriants, the value of the entries are not
          important; only the nonzero structure is used at this
          point. However for linear constraints, both the structure
          and values are important. The value of the linear
          constraints are set here and remain fixed for the remainder
          of the optimization.

          lower -> value, iteratable: The lower bounds for the constraints
          upper -> value, iteratable: The upper bounds for the constraints

        Documentation last updated:  Feb. 07, 2011 - Peter W. Jansen
        '''

        # If this is the first constraint, finalize the variables to
        # ensure no more variables can be added. 
        if self.ableToAddVariables:
            self._finalizeDesignVariables()

        if name in self.conGroupNames:
            print 'Error: The supplied name \'%s\' for a constraint group \
has already been used.'%(name)
            return
        # end if
        self.conGroupNames.add(name)

        # ------ Process the lower bound argument
        lower = numpy.atleast_1d(kwargs.pop('lower', -inf*numpy.ones(ncons)))
        if len(lower) == 1:
            lower = lower[0]*numpy.ones(ncons)
        elif len(lower) == ncons:
            pass
        else:
            print 'Error: The length of the \'lower\' argument to \
 addVarGroup is %d, but the number of variables in ncons is %d.'%(len(lower), ncons)
            sys.exit(1)
        # end if

        # ------ Process the upper bound argument
        upper = numpy.atleast_1d(kwargs.pop('upper', -inf*numpy.ones(ncons)))
        if len(upper) == 1:
            upper = upper[0]*numpy.ones(ncons)
        elif len(upper) == ncons:
            pass
        else:
            print 'Error: The length of the \'upper\' argument to \
 addVarGroup is %d, but the number of variables in ncons is %d.'%(len(upper), ncons)
            sys.exit(1)
        # end if

        # We need to do some checking if the constraint is sparse:
        if not dense:
            # First check that 'wrt' and 'jac' is given:
            if wrt is None or jac is None:
                print 'Error: \'wrt\' and \'jac\' MUST be given if the \
constriant group is sparse.'
                sys.exit(1)
            # end if

            # A special check for 'wrt' it can be 'all'; a shortcut
            # for all the variables:
            try:
                tmp = wrt.lower()
                if wrt.lower() == 'all':
                    wrt = self.variables.keys()
                else:
                    print 'Error: \'wrt\' must be a iterable list.'
                    sys.exit(1)
                # end if
            except:
                try:
                    wrt = list(wrt)
                except:
                    print 'Error: \'wrt\' must be a iterable list'
                # end try
            # end try
                    
            # Check that jac is a scipy sparse matrix:
            if sparse.issparse(jac):
                pass
            else:
                try:
                    jac = sparse.csr_matrix(jac)
                except:
                    print 'Erorr: \'jac\' is not a sparse matrix nor \
can be converted to one.'
                    sys.exit(1)
                # end try
            # end if
                    
            # Now we have to make sure the Jacobian is *actually* the
            # right size; we have supplied ncon and the wrt varibles
            # which define the size of the jacobian *should* be. 

            if ncons <> jac.shape[0]:
                print 'Error: The number of rows in \'jac\' does not equal \'ncons\''
                sys.exit(0)

            ndvs = 0
            for dvSet in wrt:
                if dvSet not in self.variables.keys():
                    print 'Error: The supplied dvSet name \'%s\' in \'wrt\' \
was not found.'%(dvSet)
                    sys.exit(1)
                else:
                    ndvs += (self.dvOffset[dvSet]['n'][1] - self.dvOffset[dvSet]['n'][0])
                # end if
            # end for

            if ndvs <> jac.shape[1]:
                print 'Error: The number of columns in \'jac\' does not equal \
the number variables defined by the sets in dvSet = %s'%(wrt)
                sys.exit(0)
            # end if
        else:
            # Create dummy values for this jacobian --- assumed
            # nonlinear so the values do not matter. 
            jac = sparse.csr_matrix(numpy.ones((ncons, self.ndvs)))
            wrt = self.variables.keys()
        # end if

        # Ensure the jacobian is CSR since the code below uses that
        # fact.
        jac = jac.tocsr()

        # For every non zero entry in jac, (ie an array of
        # len(jac.data)) then we want to know the COLUMN in the full
        # jacobain it corresonds to.
        
        jacColIndex = numpy.zeros_like(jac.indices)
        dvOffset = 0
        for dvSet in wrt:
            dvLow = self.dvOffset[dvSet]['n'][0] 
            dvHigh= self.dvOffset[dvSet]['n'][1]
            dvRange = dvHigh - dvLow

            # Loop over rows in constraint jacobian:
            for iRow in xrange(ncons):
                # Loop over the number of nonzero entries in this row:
                for ii in xrange(jac.indptr[iRow], jac.indptr[iRow+1]):
                    # ii is the index into the indices/data arrays
                    column = jac.indices[ii]
                    
                    # We only examine dvOffet + dvRange of the
                    # constriant jacobian, the remainder will be in a
                    # separate dvSet
                    if column >= dvOffset and column < dvOffset + dvRange:
                        jacColIndex[ii] = dvLow + column - dvOffset
                    # end if
                # end for
            # end for

            dvOffset += (dvHigh - dvLow)
        # end for

        self.constraints[name] = Constraint(name, dense, linear, wrt, jac, 
                                            jacColIndex, lower, upper)

        return 
        
    def reorderConstraintJacobian(self, reorder=['nonLinear','linear']):
        '''
        Here we possibly reorder the constriants to put the nonLinear
        ones first, the linear ones first or the keep the natural
        order the constraints were added. This function MUST be called
        regarless, since the con.rs and con.re values are
        computed in the function. 
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

            for iCon in self.constraints:
                con = self.constraints[iCon]
                if con.linear:
                    con.rs = rowCounter
                    rowCounter += con.ncon
                    con.re = rowCounter
        elif reorder == ['linear','nonLinear']:
            for iCon in self.constraints:
                con = self.constraints[iCon]
                if con.linear:
                    con.rs = rowCounter
                    rowCounter += con.ncon
                    con.re = rowCounter

            for iCon in self.constraints:
                con = self.constraints[iCon]
                if not con.linear:
                    con.rs = rowCounter
                    rowCounter += con.ncon
                    con.re = rowCounter
        else: # No re-ordering
            for iCon in self.constraints:
                con = self.constraints[iCon]
                con.rs = rowCounter
                rowCounter += con.ncon
                con.re = rowCounter 
            # end for
        # end if

        scale = []
        for dvSet in self.variables.keys():
            for dvGroup in self.variables[dvSet]:
                for var in self.variables[dvSet][dvGroup]:
                    scale.append(var.scale)
        self.scale = numpy.array(scale)

        return

    def processDerivatives(self, gobj_in, gcon_in, linearConstraints=False, 
                           nonlinearConstraints=True):
        '''
        This generic function is used to assemble the objective
        gradient and the constraint jacobian. The two input flags are
        used to determine which if linear/nonlinear or both are
        included. Note that all cases the size of the jacobian is
        still (ncon x ndvs), ie the full size. However, only the
        requested entries (linear/nonlinear) are included. Also note
        that this function performs the pyOpt controlled scaling that
        is transparent to the user. 
        '''
        
        # Assume gobj is already ok

        gobj = numpy.atleast_2d(gobj_in).copy()
        gobj /= self.scale

        # Data for storing the values in COOrdinate format
        data = []
        row  = []
        col  = []

        # There is also a possibility that the user has actually
        # provided jacobian that is *exactly* the correct size and
        # sparsity pattern. This will happen with pyOptSparse is used
        # in "dense" mode, ie with other pyOpt like problems. 

        allDense = True
        for iCon in self.constraints:
            if not self.constraints[iCon].dense:
                allDense = False
            # end if
        # end ofr

        # If we have all dense constraints, AND gcon_in is an array
        # AND it is the right size, we will interpret this is the
        # actual constriant jacobian and use. 
        if allDense and isinstance(gcon_in, numpy.ndarray):
            if gcon_in.shape == (self.nCon, self.ndvs):
                gcon_in[numpy.where(gcon_in==0)] = 1e-50
                # Don't forget to scale:
                tmp = gcon_in.copy()
                for i in xrange(self.ndvs):
                    tmp[:,i] /= self.scale[i]

                gcon = sparse.coo_matrix(tmp)

                return gobj, gcon
            # end if
        # end if

        # Otherwise, process constraints in the dictionary form. 
        # Loop over all constraints:
        for iCon in self.constraints:
            con = self.constraints[iCon]
            if ((linearConstraints and con.linear) or
                (nonlinearConstraints and not con.linear)):
                if not con.name in gcon_in:
                    print 'Error: the jacobian for the constraint \'%s\' was \
    not found in the returned dictionary.'%con.name
                    sys.exit(1)
                # end if
                if con.dense:
                    tmp = numpy.atleast_2d(gcon_in[iCon].copy())
                    tmp[numpy.where(tmp==0)] = 1e-50
                    tmp = sparse.csr_matrix(tmp)
                else:
                    tmp = sparse.csr_matrix(gcon_in[iCon].copy())
                # end if

                if tmp.shape <> con.jac.shape:
                    print 'Error: The jacobian for constraint group \'%s\' \
was not the correct size. The supplied jacobian has shape %s, but must be \
shape %s.'%(con.name, gcon_in[iCon].shape, con.jac.shape)
                    sys.exit(1)
                # end if

                if tmp.nnz <> con.jac.nnz:
                    print 'Error: The number of nonzero elements for \
 constraint group \'%s\' was not the correct size. The supplied jacobian has \
%d nonzero entries, but must contain %d nonzero entries.'%(
                        con.name, tmp.nnz, con.jac.nnz)
                    sys.exit(1)
                # end if
                
                # Loop over rows in constraint jacobian:
                for iRow in xrange(con.ncon):
                    # Loop over the number of nonzero entries in this row:
                    for ii in xrange(con.jac.indptr[iRow], con.jac.indptr[iRow+1]):
                        row.append(con.rs + iRow)
                        col.append(con.jacColIndex[ii])
                        data.append(tmp.data[ii]/self.scale[con.jacColIndex[ii]])
                    # end for
                # end for
            # end if (nonlinear)
        # end for

        gcon = sparse.coo_matrix((data, (row, col)),(self.nCon, self.ndvs))

        return gobj, gcon

    def processConstraints(self, tmp, linearConstraints=False, 
                           nonlinearConstraints=True):
        '''
        Assemble the constraint vector from the returned dictionary
        '''

        # We will actually be a little leniant here; the user CAN
        # return an iterable of the correct length and we will accept
        # that. Otherwise we will use the dictionary formulation
        error = False

        if not isinstance(tmp, dict):
            fcon = numpy.atleast_1d(tmp)
            if len(fcon) == self.nnCon:
                return fcon
            else:
                print 'Error: The constraint array was the incorrect size. \
It must contain %d elements (nonlinear constraints only), but an arrary of \
size %d was given.'%(self.nnCon, len(fcon))
                error = True
            # end if
        else:
            # Process as a dictionary:
            # Loop over (nonlinear) constraints and extract as required:
            fcon = []
            for iCon in self.constraints:
                if not self.constraints[iCon].linear:
                    if iCon in tmp:
                        # Make sure it is at least 1dimension:
                        c = numpy.atleast_1d(tmp[iCon])
                        
                        # Make sure it is the correct size:
                        if len(c) == self.constraints[iCon].ncon:
                            fcon.extend(c)
                        else:
                            print 'Error: %d constraint values were returned \
    %s, but expected %d.'%(len(tmp[iCon]), iCon, self.variables[iCon].ncon)
                            error = True
                        # end if
                    else:
                        print 'Error: No constraint values were found for the \
constraint %s.'%(iCon)
                        error = True
                    # end if
                # end if
            # end for
            # Finally convert to array:
            fcon = numpy.array(fcon)
        # end if
 
        if error:
            sys.exit(1)

        return fcon

    def convertToDense(self):
        '''
        Take a sparse optimization problem definition and convert to a
        dense representation for use in the rest of pyOpt
        '''

        # Variables are the same except the stupid underscore
        self._variables = {}
        ii = 0
        for dvSet in self.variables.keys():
            for dvGroup in self.variables[dvSet]:
                for i in xrange(len(self.variables[dvSet][dvGroup])):
                    self._variables[ii] = self.variables[dvSet][dvGroup][i]
                    ii += 1
                # end for
            # end for
        # end for
        
        from pyOpt import Constraint as pyOptConstraint

        # Constraints have to be done individually
        self._constraints = {}
        ii = 0
        for iCon in self.constraints:
            con = self.constraints[iCon]
            for i in xrange(con.ncon):
                self._constraints[ii] = pyOptConstraint(
                    con.name+'_%d'%(i), type=con.type, lower=con.lower[i],
                    upper=con.upper[i])
                ii += 1
            # end for
        # end for
        
        self._objectives = {}
        ii = 0
        for obj in self.objectives:
            self._objectives[ii] = self.objectives[obj]

        self.assembleFullConstraintJacobian(reorder=None)

        return

    def processX(self, x):
        '''
        Take the flattened array of design variables and return a dict
        '''
        if self.use_groups:
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

    def __str__(self):
        
        '''
        Print Structured Optimization Problem
        
        Documentation last updated:  April. 30, 2008 - Peter W. Jansen
        '''
        
        text = '''\nOptimization Problem -- %s\n%s\n
        Objective Function: %s\n\n    Objectives:
        Name        Value        Optimum\n''' %(self.name,'='*80,self.obj_fun.__name__)

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
        sparsity structure. This helps the user visualize what pyOpt
        has been given and helps ensure it is what the user expected. 
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
            maxConNameLen = max(maxConNameLen, len(self.constraints[iCon].name))
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
            varCenters.append(nCol + len(iVar)/2 + 1)
            nCol += len(iVar) 
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
            l = len(iVar)
            txt[0, iCol+1 :iCol + l+1] = list(iVar)
            txt[2:-1, iCol + l + 2] = '|'
            iCol += l + 3

        # Print the constraint names;
        iRow = 2

        # Do the nonlinear ones first:
        for iCon in self.constraints:
            if not self.constraints[iCon].linear:
                name = self.constraints[iCon].name
                l = len(name)
                # The name
                txt[iRow, maxConNameLen-l:maxConNameLen] = list(name)

                # Now we write a 'D' for dense, 'S' for sparse or nothing. 
                for iVar in xrange(len(self.variables)):

                    if self.constraints[iCon].dense:
                        txt[iRow, varCenters[iVar]] = 'D'
                    elif self.variables.keys()[iVar] in self.constraints[iCon].wrt:
                        txt[iRow, varCenters[iVar]] = 'S'
                    # end if
                # end for

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
            if self.constraints[iCon].linear:
                name = self.constraints[iCon].name
                l = len(name)
                # The name
                txt[iRow, maxConNameLen-l:maxConNameLen] = list(name)

                # Now we write a 'D' for dense, 'S' for sparse or nothing. 
                for iVar in xrange(len(self.variables)):

                    if self.constraints[iCon].dense:
                        txt[iRow, varCenters[iVar]] = 'D'
                    elif self.variables.keys()[iVar] in self.constraints[iCon].wrt:
                        txt[iRow, varCenters[iVar]] = 'S'
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
    
    
