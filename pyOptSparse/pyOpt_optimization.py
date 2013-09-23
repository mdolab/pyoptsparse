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
    - add variable group error when groups have the same name
    - add method for addVar2Group
    - pickle wrapping ?!
    - save class __str__ info to file (text/TeX) ?
    - warm start from other opts?
    - class for core sensitivity?
    - class for history?
'''

# =============================================================================
# Standard Python modules
# =============================================================================
import os, sys
from collections import OrderedDict

# =============================================================================
# External Python modules
# =============================================================================
import numpy
from scipy import sparse
# =============================================================================
# Extension modules
# =============================================================================
from pyOpt import Variable
from pyOpt import Objective
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


        # A set to keep track of user-supplied names -- don't let the
        # user reuse names for design varible sets, design variable
        # groups, or constraint groups
        self.allNames = set()

        # Flag to determine if adding variables is legal. 
        self.ableToAddVariables = True

        # Add a default variable set:
        #Qself.addVarSet('default')

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
                self.dvOffset[dvSet][dvGroup] = [dvCounter, dvCounter + n]
                dvCounter += n
            self.dvOffset[dvSet]['n'][1] = dvCounter
        # end for
        self.ndvs = dvCounter

        return

    def addVarSet(self, name):
        '''An outer grouping of design variables. These sets are used
        when specifiying the sparsity structuer of the constraint
        jacobian
        '''

        self._checkOkToAddVariables()
        
        if name in self.allNames:
            print 'Error: The supplied name \'%s\' for a variable set \
has already been used.'%(name)
            return
        # end if
        self.allNames.add(name)
        self.variables[name]=OrderedDict()

        return

    def addVar(self, name, *args, **kwargs):
        '''
        convenience function. See addVarGroup for more information
        '''

        self.addVarGroup(name, 1, *args, **kwargs)

        return 

    def addVarGroup(self, name, nvars, type='c', value=0.0, varSet='default', **kwargs):
        
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

        if name in self.allNames:
            print 'Error: The supplied name \'%s\' for a variable group \
has already been used.'%(name)
            return
        # end if
        self.allNames.add(name)

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

        # Now create all the variable objects
        self.variables[varSet][name] = []
        for iVar in xrange(nvars):
            varName = name + '_%d'%(iVar)
            self.variables[varSet][name].append(
                Variable(varName, type=type, value=value[iVar], 
                         lower=lower[iVar], upper=upper[iVar]))
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
          wrt to variables 'A' and the last 5 columes wrt variables
          'B'. If wrt was given as wrt=['B','A'] the opposite order
          would be used. Note that the wrt order MUST be used in all
          subequent returns of the jacobian. 

          jac -> scipy.sparse matrix: The sparse matrix jacobian. For
          nonlinear constriants, the value of the entries are not
          important; only the structure is used at this point. However
          for linear constraints, both the structure and values are
          important. The value of the linear constraints are set here
          and remain fixed for the remainder of the optimization. 

          lower -> value, iteratable: The lower bounds for the constraints
          upper -> value, iteratable: The upper bounds for the constraints
        
        Documentation last updated:  Feb. 07, 2011 - Peter W. Jansen
        '''

        if self.ableToAddVariables:
            self._finalizeDesignVariables()

        if name in self.allNames:
            print 'Error: The supplied name \'%s\' for a constraint group \
has already been used.'%(name)
            return
        # end if
        self.allNames.add(name)

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
                    jac = sparse.lil_matrix(jac)
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

            # We now know the jac is the right size we now want to
            # section it according to the dvSets given

            iStart = 0
            cs = [] # Column start in full jacobian
            ce = [] # Column end in full jacobian
            jcs = [] # Column start in this constraint jacobian
            jce = [] # Column end in this constraint jacobian
            for dvSet in wrt:
                cs.append(self.dvOffset[dvSet]['n'][0] )
                ce.append(self.dvOffset[dvSet]['n'][1])

                jcs.append(iStart)
                iEnd = iStart + (ce[-1] - cs[-1])
                jce.append(iEnd)

                iStart = iEnd
            # end if
        else:
            # Create dummy values for this jacobian --- assumed
            # nonlinear so the values do not matter. 
            jac = numpy.ones((ncons, self.ndvs))
            cs = [0]
            ce = [self.ndvs]
            jcs = [0]
            jce = [self.ndvs]
        # end if

        self.constraints[name] = Constraint(name, dense, linear, wrt, jac, 
                                            cs, ce, jcs, jce, lower, upper)

        return 
        
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
    
    
