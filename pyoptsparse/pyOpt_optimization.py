#!/usr/bin/env python
"""
pyOptSparse_optimization

Holds the Python Design Optimization Class

The main purpose, of this class is to describe the struture and
potentially, sparsity pattern of an optimization problem.

Copyright (c) 2013-2014 by Dr. Gaetan Kenway
All rights reserved.

Developers:
-----------
- Dr. Gaetan K.W. Kenway (GKK)

History
-------
    v. 1.0  - Initial Class Creation (GKK, 2013)
"""
from __future__ import print_function

# =============================================================================
# Standard Python modules
# =============================================================================
import copy
import os
import shelve
try:
    from collections import OrderedDict
except ImportError:
    try:
        from ordereddict import OrderedDict
    except ImportError:
        print('Could not find any OrderedDict class. For 2.6 and earlier, \
use:\n pip install ordereddict')
    
# =============================================================================
# External Python modules
# =============================================================================
import numpy

# =============================================================================
# Extension modules
# =============================================================================
from .pyOpt_variable import Variable
from .pyOpt_objective import Objective
from .pyOpt_constraint import Constraint
from .pyOpt_error import Error
from .pyOpt_MPI import MPI
from .pyOpt_utils import *
# =============================================================================
# Misc Definitions
# =============================================================================
INFINITY = 1e20 
# =============================================================================
# Optimization Class
# =============================================================================
class Optimization(object):
    """
    Create a description of an optimization probelem. 

    Parameters
    ----------
    name : str
        Name given to optimization problem. This is name is currently
        not used for anything, but may be in the future.

    objFun : python function
        Python function handle of function used to evaluate the objective
        function.

    comm : MPI intra communication
        The communicator this problem will be solved on. This is
        required for both analysis when the objective is computed in
        parallel as well as to use the internal parallel gradient
        computations. Defaults to MPI.COMM_WORLD if not given. 
        """
      
    def __init__(self, name, objFun, comm=None):

        self.name = name
        self.objFun = objFun
        if comm is None:
            self.comm = MPI.COMM_WORLD
        else:
            self.comm = comm
            
        # Ordered dictionaries to keep track of variables and constraints
        self.variables = OrderedDict()
        self.constraints = OrderedDict()
        self.objectives = OrderedDict()
        self.dvOffset =  OrderedDict()

        # Flag to determine if adding variables is legal. 
        self.denseJacobianOK = True
        
        # Variables to be set in finalizeConstraints
        # have finalized the specification of the variable and the
        # constraints
        self.ndvs = None
        self.conScale = None
        self.nCon = None
        self.xscale = None
        self.linearJacobian = None
        self.dummyConstraint = False
        
    def addVar(self, name, *args, **kwargs):
        """
        This is a convience function. It simply calls addVarGroup()
        with nVars=1. Variables added with addVar() are returned as
        *scalars*. 
        """
        self.addVarGroup(name, 1, *args, scalar=True, **kwargs)

    def checkVarName(self, varName):
        """
        Check if the desired variable name varName if has already been
        added. If it is has already been added, return a managled name
        (a number appended) that *is* valid. This is intended to be used
        by classes that automatically add variables to pyOptSparse

        Parameters
        ----------
        varName : str
            Variable name to check validity on

        Returns
        -------
        validName : str
            A valid variable name. May be the same as varName it that
            was, in fact, a valid name.
            """
        if varName not in self.variables:
            return varName
        else:
            i = 0
            validName = varName + '_%d'% i
            while validName in self.variables:
                i += 1
                validName = varName + '_%d'% i
            return validName

    def checkConName(self, conName):
        """
        Check if the desired constraint name has already been
        added. If it is has already been added, return a managled name
        (a number appended) that *is* valid. This is intended to be used
        by classes that automatically add constraints to pyOptSparse.

        Parameters
        ----------
        conName : str
           Constraint name to check validity on

        Returns
        -------
        validName : str
            A valid constraint name. May be the same as conName it that
            was, in fact, a valid name.
            """
        if conName not in self.constraints:
            return conName
        else:
            i = 0
            validName = conName + '_%d'% i
            while validName in self.constraints:
                i += 1
                validName = conName + '_%d'% i
            return validName

    def addVarGroup(self, name, nVars, type='c', value=0.0, 
                    lower=None, upper=None, scale=1.0, 
                    choices=None, **kwargs):
        """
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
            
        scale : scalar or array.  Define a user supplied scaling
            variable for the design variable group.  This is often
            necessary when design variables of widely varraying
            magnitudes are used within the same
            optimization. Scalar/array usage is the same as value
            keyword.

        choices : list
            Specify a list of choices for discrete design variables
            
        Examples
        --------
        >>> # Add a single design variable 'alpha'
        >>> optProb.addVar('alpha', type='c', value=2.0, lower=0.0, upper=10.0, \
        scale=0.1)
        >>> # Add 10 unscaled variables of 0.5 between 0 and 1 with name 'y'
        >>> optProb.addVarGroup('y', type='c', value=0.5, lower=0.0, upper=1.0, \
        scale=1.0)
        
        Notes
        -----
        Calling addVar() and addVarGroup(..., nVars=1, ...) are
        **NOT** equilivant! The variable added with addVar() will be
        returned as scalar, while variable returned from addVarGroup
        will be an array of length 1.

        It is recommended that the addVar() and addVarGroup() calls
        follow the examples above by including all the keyword
        arguments. This make it very clear the itent of the script's
        author. The type, value, lower, upper and scale should be
        given for all variables even if the default value is used. 
        """

        # Check that the nVars is > 0. 
        if nVars < 1:
            raise Error("The 'nVars' argument to addVarGroup must be greater "
                        "than or equal to 1. The bad DV is %s."%name)

        # Check that the type is ok:
        if type not in ['c', 'i', 'd']:
            raise Error("Type must be one of 'c' for continuous, "
                        "'i' for integer or 'd' for discrete.")
                    
        # ------ Process the value arguement
        value = numpy.atleast_1d(value).real
        if len(value) == 1:
            value = value[0]*numpy.ones(nVars)
        elif len(value) == nVars:
            pass
        else:
            raise Error("The length of the 'value' argument to "
                        "addVarGroup is %d, but the number of "
                        "variables in nVars is %d."% (len(value), nVars))

        if lower is None:
            lower = [None for i in range(nVars)]
        elif numpy.isscalar(lower):
            lower = lower*numpy.ones(nVars)
        elif len(lower) == nVars:
            lower = numpy.atleast_1d(lower).real
        else:
            raise Error("The 'lower' argument to addVarGroup is "
                        "invalid. It must be None, a scalar, or a "
                        "list/array or length nVars=%d." %(nVars))

        if upper is None:
            upper = [None for i in range(nVars)]
        elif numpy.isscalar(upper):
            upper = upper*numpy.ones(nVars)
        elif len(upper) == nVars:
            upper = numpy.atleast_1d(upper).real
        else:
            raise Error("The 'upper' argument to addVarGroup is "
                        "invalid. It must be None, a scalar, or a "
                        "list/array or length nVars=%d." %(nVars))

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
                raise Error("The length of the 'scale' argument to "
                            "addVarGroup is %d, but the number of "
                            "variables in nVars is %d."% (
                                len(scale), nVars))

        # Determine if scalar i.e. it was called from addVar():
        scalar = kwargs.pop('scalar', False)

        # Now create all the variable objects
        varList = []
        for iVar in range(nVars):
            varName = name + '_%d'% iVar
            varList.append(Variable(varName, type=type, value=value[iVar],
                                    lower=lower[iVar], upper=upper[iVar],
                                    scale=scale[iVar], scalar=scalar, 
                                    choices=choices))

        if name in self.variables:
            # Check that the variables happen to be the same
            err = False
            if not len(self.variables[name]) == len(varList):
                raise Error("The supplied name '%s' for a variable group "
                            "has already been used!"% name)
            for i in range(len(varList)):
                if not varList[i] == self.variables[name][i]:
                    raise Error("The supplied name '%s' for a variable group "
                                "has already been used!"% name)
            # We we got here, we know that the variables we wanted to
            # add are **EXACTLY** the same so that's cool. We'll just
            # overwrite with the varList below. 
        else:
            # Finally we set the variable list
            self.variables[name] = varList

    def delVar(self, name):
        """
        Delete a variable or variable group

        Parameters
        ----------
        name : str
           Name of variable or variable group to remove
           """
        try:
            self.variables.pop(name)
        except KeyError:
            print('%s was not a valid design variable name.'% name)
            
    def _reduceDict(self, variables):
        """
        This is a specialized function that is used to communicate
        variables from dictionaries across the comm to ensure that all
        processors end up with the same dictionary. It is used for
        communicating the design variables and constraints, which may
        be specified on different processors independently.
        """
        
        # Step 1: Gather just the key names:
        allKeys = self.comm.gather(list(variables.keys()), root=0)
       
        # Step 2: Determine the unique set:
        procKeys = {}
        if self.comm.rank == 0:
            # We can do the reduction efficiently using a dictionary: The
            # algorithm is as follows: . Loop over the processors in order,
            # and check if key is in procKeys. If it isn't, add with proc
            # ID. This ensures that when we're done, the keys of 'procKeys'
            # contains all the unique values we need, AND it has a single
            # (lowest proc) that contains that key
            for iProc in range(len(allKeys)):
                for key in allKeys[iProc]:
                    if not key in procKeys:
                        procKeys[key] = iProc

            # Now pop any keys out with iProc = 0, since we want the
            # list of ones NOT one the root proc
            for key in list(procKeys.keys()):
                if procKeys[key] == 0:
                    procKeys.pop(key)

        # Step 3. Now broadcast this back to everyone
        procKeys = self.comm.bcast(procKeys, root=0)
     
        # Step 4. The required processors can send the variables
        if self.comm.rank == 0:
            for key in procKeys:
                variables[key] = self.comm.recv(source=procKeys[key], tag=0)
        else:
            for key in procKeys:
                if procKeys[key] == self.comm.rank:
                    self.comm.send(variables[key], dest=0, tag=0)

        # Step 5. And we finally broadcast the final list back:
        variables = self.comm.bcast(variables, root=0)

        return variables

    def addObj(self, name, *args, **kwargs):
        """
        Add Objective into Objectives Set
        """
        
        self.objectives[name] = Objective(name, *args, **kwargs)

    def addCon(self, name, *args, **kwargs):
        """
        Convenience function. See addConGroup() for more information
        """
        
        self.addConGroup(name, 1, *args, **kwargs)

    def addConGroup(self, name, nCon, lower=None, upper=None, scale=1.0, 
                    linear=False, wrt=None, jac=None):
        """Add a group of variables into a variable set. This is the main
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
            specifies for what dvs have non-zero jacobian values
            for this set of constraints. The order is not important.

        jac : dictionary
            For linear and sparse non-linear constraints, the constraint
            jacobian must be passed in. The structure is jac dictionary
            is as follows:

            {'dvName1':<matrix1>, 'dvName2', <matrix1>, ...}

            They keys of the jacobian must correpsond to the dvGroups
            given in the wrt keyword argument. The dimensions of each
            "chunk" of the constraint jacobian must be consistent. For
            example, <matrix1> must have a shape of (nCon, nDvs) where
            nDVs is the total number of design variables in
            dvName1. <matrix1> may be a dense numpy array or it may be
            scipy sparse matrix. However, it is *HIGHLY* recommended
            that sparse constraints are supplied to pyOptSparse using
            the pyOptSparse's simplified sparse matrix format. The
            reason for this is that it is *impossible* for force scipy
            sparse matrices to keep a fixed sparsity pattern; if the
            sparsity pattern changes during an optimization, *IT WILL
            FAIL*.

            The three simplified pyOptSparse sparse matrix formats are
            sumarized below:

            mat = {'coo':[row, col, data], 'shape':[nrow, ncols]} # A coo matrix
            mat = {'csr':[rowp, colind, data], 'shape':[nrow, ncols]} # A csr matrix
            mat = {'coo':[colp, rowind, data], 'shape':[nrow, ncols]} # A csc matrix

            Note that for nonlinear constraints (linear=False), the
            values themselves in the matrices in jac do not matter, 
            but the sparsity structure **does** matter. It is
            imparative that entries that will at some point have
            non-zero entries have non-zero entries in jac
            argument. That is, we do not let the sparsity structure of
            the jacobian change throughout the optimization. This
            stipulation is automatically checked internally.

        """

        if name in self.constraints:
            raise Error("The supplied name '%s' for a constraint group "
                        "has already been used."% name)

        # Simply add constraint object
        self.constraints[name] = Constraint(
            name, nCon, linear, wrt, jac, lower, upper, scale)

    def getDVs(self):
        """
        Return a dictionary of the design variables. In most common
        usage, this function is not required.

        Returns
        -------
        outDVs : dict
            The dictionary of variables. This is the same as 'x' that
            would be used to call the user objective function. 
        """ 

        outDVs = {}
        for dvGroup in self.variables:
            nvar = len(self.variables[dvGroup])
            # If it is a single DV, return a scalar rather than a numpy array
            if nvar == 1:
                var = self.variables[dvGroup][0]
                outDVs[dvGroup] = var.value/var.scale
            else:
                outDVs[dvGroup] = numpy.zeros(nvar)
                for i in range(nvar):
                    var = self.variables[dvGroup][i]
                    outDVs[dvGroup][i] = var.value/var.scale

        return outDVs

    def setDVs(self, inDVs):
        """
        set the problem design variables from a dictionary. In most
        common usage, this function is not required.
        
        Parameters
        ----------
        inDVs : dict
            The dictionary of variables. This dictionary is like the
            'x' that would be used to call the user objective
            function.
        """
        for dvGroup in self.variables:
            if dvGroup in inDVs:
                nvar = len(self.variables[dvGroup])
                for i in range(nvar):
                    var = self.variables[dvGroup][i]
                    if numpy.isscalar(inDVs[dvGroup]):
                        var.value = inDVs[dvGroup]*var.scale
                    else:
                        var.value = inDVs[dvGroup][i]*var.scale

    def setDVsFromHistory(self, histFile, key=None):
        """
        Set optimization variables from a previous optimization. This
        is like a cold start, but some variables may have been added
        or removed from the previous optimization. This will try to
        set all variables it can.

        Parameters
        ----------
        histFile : str
            Filename of the history file to read
        key : str
            Key of the history file to use for the x values. The
            default is None which will use the last x-value stored in
            the dictionary.
        """
        
        if os.path.exists(histFile):
            hist = shelve.open(histFile, flag='r')
            if key is None:
                key = hist['last']
                
            self.setDVs(hist[key]['xuser'])
            hist.close()
        else:
            raise Error("History file '%s' not found!."% histFile)
        
    def printSparsity(self):
        """
        This function prints an (ascii) visualization of the jacobian
        sparsity structure. This helps the user visualize what
        pyOptSparse has been given and helps ensure it is what the
        user expected. It is highly recommended this function be
        called before the start of every optimization to verify the
        optimization problem setup.

        Warnings
        --------
        This function is **collective** on the optProb comm. It is
        therefore necessary to call this function on **all**
        processors of the optProb comm. 

        """
        self.finalizeDesignVariables()
        self.finalizeConstraints()

        if self.comm.rank != 0:
            return
    
        # Header describing what we are printing:
        print('+'+'-'*78+'-'+'+')
        print('|' + ' '*19 +'Sparsity structure of constraint Jacobian' + ' '*19 + '|')
        print('+'+'-'*78+'-'+'+')

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
            maxConNameLen = max(maxConNameLen,
                                len(con.name)+6+int(numpy.log10(con.ncon))+1)
            nRow += 1 # Line

        # And now the columns:
        nCol = maxConNameLen
        nCol += 2 # Space plus line
        varCenters = []
        for dvGroup in self.variables:
            nvar = self.dvOffset[dvGroup][1] - self.dvOffset[dvGroup][0]
            var_str = dvGroup + ' (%d)'% nvar

            varCenters.append(nCol + len(var_str)/2 + 1)
            nCol += len(var_str)
            nCol += 2 # Spaces on either side
            nCol += 1 # Line 

        txt = numpy.zeros((nRow, nCol), dtype=str)
        txt[:, :] = ' '
        # Outline of the matrix on left and top
        txt[1, maxConNameLen+1:-1] = '-'
        txt[2:-1, maxConNameLen+1] = '|'
     
        # Print the variable names:
        iCol = maxConNameLen + 2
        for dvGroup in self.variables:
            nvar = self.dvOffset[dvGroup][1] - self.dvOffset[dvGroup][0]
            var_str = dvGroup + ' (%d)'% nvar
            l = len(var_str)
            txt[0, iCol+1 :iCol + l+1] = list(var_str)
            txt[2:-1, iCol + l + 2] = '|'
            iCol += l + 3

        # Print the constraint names;
        iRow = 2

        for iCon in self.constraints:
            con = self.constraints[iCon]
            name = con.name
            if con.linear:
                name = name + '(L)'

            name = name + ' (%d)'% con.ncon
            l = len(name)
            # The name
            txt[iRow, maxConNameLen-l:maxConNameLen] = list(name)

            # Now we write a 'X' if there is something there:
            varKeys = list(self.variables.keys())
            for dvGroup in range(len(varKeys)):
                if varKeys[dvGroup] in con.wrt:
                    txt[iRow, varCenters[dvGroup]] = 'X'

            # The separator
            txt[iRow+1, maxConNameLen+1:] = '-'
            iRow += 2
            
        # Corners - just to make it nice :-)
        txt[1, maxConNameLen+1] = '+'
        txt[-1, maxConNameLen+1] = '+'
        txt[1, -1] = '+'
        txt[-1, -1] = '+'

        for i in range(len(txt)):
            print(''.join(txt[i]))

    def getDVConIndex(self, startIndex=1, printIndex=True):
        """
        Return the index of a scalar DV/constraint, or the beginning
        and end index (inclusive) of a DV/constraint array.
        This is useful for looking at SNOPT gradient check output, 
        and the default startIndex=1 is for that purpose
        """

        # Get the begin and end index (inclusive) of design variables
        # using infomation from finalizeDesignVariables()
        dvIndex = OrderedDict()
        # Loop over the actual DV names 
        for dvGroup in self.dvOffset:
            ind0 = self.dvOffset[dvGroup][0] + startIndex
            ind1 = self.dvOffset[dvGroup][1] + startIndex
            # if it is a scalar DV, return just the index
            if ind1 - ind0 == 1:
                dvIndex[dvGroup] = [ind0]
            else:
                dvIndex[dvGroup] = [ind0, ind1-1]

        # Get the begin and end index (inclusive) of constraints
        conIndex = OrderedDict()
        conCounter = startIndex
        for iCon in self.constraints:
            n = self.constraints[iCon].ncon
            if n == 1:
                conIndex[iCon] = [conCounter]
            else:
                conIndex[iCon] = [conCounter, conCounter+n-1]
            conCounter += n

        # Print them all to terminal
        if printIndex and self.comm.rank == 0:
            print('### DESIGN VARIABLES ###')
            for dvGroup in dvIndex:
                print(dvGroup, dvIndex[dvGroup])
            print('### CONSTRAINTS ###')
            for conKey in conIndex:
                print(conKey, conIndex[conKey])
        
        return dvIndex, conIndex

#=======================================================================
#       All the functions from here down should not need to be called
#       by the user. Most functions are public since the individual
#       optimizers need to be able to call them
#=======================================================================

    def finalizeDesignVariables(self):
        """
        Communicate design variables potentially from different
        processors and form the DVOffset dict. This routine should be
        called from the individual optimizers
        """

        # First thing we need is to determine the consistent set of
        # variables from all processors.
        self.variables = self._reduceDict(self.variables)

        dvCounter = 0
        self.dvOffset = OrderedDict()

        for dvGroup in self.variables:
            n = len(self.variables[dvGroup])
            self.dvOffset[dvGroup] = [
                dvCounter, dvCounter + n, 
                self.variables[dvGroup][0].scalar]
            dvCounter += n
        self.ndvs = dvCounter
    
    def finalizeConstraints(self):
        """
        **This function should not need to be called by the user**

        There are several functions for this routine:

        1. Determine the number of constraints
        
        2. Determine the final scaling array for the design variables

        3. Determine if it is possible to return a complete dense
           jacobian. Most of this time, we should be using the dictionary-
           based return
        """

        # First thing we need is to determine the consistent set of
        # constraints from all processors
        self.constraints = self._reduceDict(self.constraints)

        # ----------------------------------------------------
        # Step 1. Determine number of constraints and scaling:
        # ----------------------------------------------------

        # Determine number of constraints
        self.nCon = 0
        for iCon in self.constraints:
            self.nCon += self.constraints[iCon].ncon
            
        # Loop over the constraints assigning the row start (rs) and
        # row end (re) values. The actual ordering depends on if
        # constraints are reordered or not.
        rowCounter = 0
        conScale = numpy.zeros(self.nCon)
        for iCon in self.constraints:
            con = self.constraints[iCon]
            con.finialize(self.variables, self.dvOffset, rowCounter)
            rowCounter += con.ncon
            conScale[con.rs:con.re] = con.scale

        if self.nCon > 0:
            self.conScale = conScale
        else:
            self.conScale = None
        
        # -----------------------------------------
        # Step 2. Assemble design variable scaling
        # -----------------------------------------
        xscale = []
        for dvGroup in self.variables:
            for var in self.variables[dvGroup]:
                xscale.append(var.scale)
        self.invXScale = 1.0/numpy.array(xscale)

        # ----------------------------------------
        # Step 3. Determine if dense return is OK
        # ----------------------------------------
        allVars = set(self.variables.keys())
        for iCon in self.constraints:
            con = self.constraints[iCon]
            # All entries of con.wrt in allVars
            if not set(con.wrt) <= allVars:
                # If any constrant 'wrt' is not fully in allVars we
                # can't do a dense return
                self.denseJacobianOK = False

        # ---------------------------------------------
        # Step 3. Final jacobian for linear constraints
        # ---------------------------------------------
        for iCon in self.constraints:
            con = self.constraints[iCon]
            if con.linear:
                data = []
                row = []
                col = []

                for dvGroup in con.jac:
                    # ss means 'start - stop'
                    ss = self.dvOffset[dvGroup]

                    row.extend(con.jac[dvGroup]['coo'][IROW])
                    col.extend(con.jac[dvGroup]['coo'][ICOL] +ss[0])
                    data.extend(con.jac[dvGroup]['coo'][IDATA])
                    
                # Now create a coo, convert to CSR and store
                con.linearJacobian = {'coo':[row, col, data],
                                      'shape':[con.ncon, self.ndvs]}

    def getOrdering(self, conOrder, oneSided, noEquality=False):
        """
        Internal function that is used to produce a index list that 
        reorders the constraints the way a particular optimizer needs. 

        Parameters
        ----------
        conOrder : list
            This must contain the following 4 strings: 'ni', 'li',
            'ne', 'le' whcih stand for nonlinear inequality, linear
            inequality, nonlinear equality and linear equality. This
            defines the order that the optimizer wants the constraints

        oneSided : bool
           Flag to do all constraints as one-sided instead of two
           sided. Most optimzers need this but some can deal with the
           two-sided constraints properly (snopt and ipopt for
           example)

        noEquality : bool
           Flag to split equality constraints into two inequality
           constraints. Some optimizers (CONMIN for example) can't do
           equality constraints explictly. 
           """

        # Now for the fun part determine what *actual* order the
        # constraints need to be in: We recongize the following
        # constraint types:
        # ne : nonlinear equality
        # ni : nonlinear inequality
        # le : linear equality
        # li : linear inequality

        # The oneSided flag determines if we use the one or two sided
        # constraints. The result of the following calculation is the
        # a single index vector that that maps the natural ordering of
        # the constraints to the order that optimizer has
        # requested. This will be returned so the optimizer can do
        # what they want with it.

        if self.nCon == 0:
            if self.dummyConstraint:
                return [], [-INFINITY], [INFINITY], None
            else:
                return numpy.array([], 'd')

        indices = []
        fact = []
        lower = []
        upper = []
        
        for conType in conOrder:
            for iCon in self.constraints:
                con = self.constraints[iCon]
                # Make the code below easier to read:
                econ = con.equalityConstraints
                if oneSided:
                    icon = con.oneSidedConstraints
                else:
                    icon = con.twoSidedConstraints
              
                if conType == 'ne' and not con.linear:
                    if noEquality:
                        # Expand Equality constraint to two:
                        indices.extend(con.rs + econ['ind'])
                        fact.extend(econ['fact'])
                        lower.extend(econ['value'])
                        upper.extend(econ['value'])
                        #....And the other side
                        indices.extend(con.rs + econ['ind'])
                        fact.extend(-1.0*econ['fact'])
                        lower.extend(econ['value'])
                        upper.extend(econ['value'])

                    else:
                        indices.extend(con.rs + econ['ind'])
                        fact.extend(econ['fact'])
                        lower.extend(econ['value'])
                        upper.extend(econ['value'])
                    
                if conType == 'ni' and not con.linear:
                    indices.extend(con.rs + icon['ind'])
                    fact.extend(icon['fact'])
                    lower.extend(icon['lower'])
                    upper.extend(icon['upper'])
                        
                if conType == 'le' and con.linear:
                    if noEquality:
                        # Expand Equality constraint to two:
                        indices.extend(con.rs + econ['ind'])
                        fact.extend(econ['fact'])
                        lower.extend([-INFINITY]*len(econ['fact']))
                        upper.extend(econ['value'])
                        #....And the other side
                        indices.extend(con.rs + econ['ind'])
                        fact.extend(-1.0*econ['fact'])
                        lower.extend([-INFINITY]*len(econ['fact']))
                        upper.extend(-econ['value'])
                    else:
                        indices.extend(con.rs + econ['ind'])
                        fact.extend(econ['fact'])
                        lower.extend(econ['value'])
                        upper.extend(econ['value'])

                if conType == 'li' and con.linear:
                    indices.extend(con.rs + icon['ind'])
                    fact.extend(icon['fact'])
                    lower.extend(icon['lower'])
                    upper.extend(icon['upper'])

        if len(fact) == 0:
            fact = None
        return numpy.array(indices), numpy.array(lower), numpy.array(upper), fact
    
    def processX(self, x):
        """
        **This function should not need to be called by the user**

        Take the flattened array of variables in 'x' and return a
        dictionary of variables keyed on the name of each variable. 

        Parameters
        ----------
        x : array
            Flattened array from optimizer
        """
        xg = {}
        imax = 0
        for dvGroup in self.variables:
            istart = self.dvOffset[dvGroup][0]
            iend = self.dvOffset[dvGroup][1]
            scalar = self.dvOffset[dvGroup][2]
            imax = max(imax, iend)
            try:
                if scalar:
                    xg[dvGroup] = x[istart]
                else:
                    xg[dvGroup] = x[istart:iend].copy()
            except IndexError:
                raise Error("Error processing x. There "
                            "is a mismatch in the number of variables.")
        if imax != self.ndvs:
            raise Error("Error processing x. There "
                        "is a mismatch in the number of variables.")
        return xg

    def deProcessX(self, x):
        """
        **This function should not need to be called by the user**

        Take the dictionary form of x and convert back to flattened
        array.

        Parameters
        ----------
        x : dict
            Dictionary form of variables

        Returns
        -------
        x_array : array
            Flattened array of variables
        """

        x_array = numpy.zeros(self.ndvs)
        imax = 0
        for dvGroup in self.variables:
            istart = self.dvOffset[dvGroup][0]
            iend = self.dvOffset[dvGroup][1]
            imax = max(imax, iend)
            scalar = self.dvOffset[dvGroup][2]
            try:
                if scalar:
                    x_array[istart] = x[dvGroup]
                else:
                    x_array[istart:iend] = x[dvGroup]
            except IndexError:
                raise Error("Error deprocessing x. There "
                            "is a mismatch in the number of variables.")
        if imax != self.ndvs:
            raise Error("Error deprocessing x. There is a mismatch in the"
                        " number of variables.")
            
        return x_array

    def processObjective(self, funcs, scaled=True):
        """
        **This function should not need to be called by the user**

        This is currently just a stub-function. It is here since it
        the future we may have to deal with multiple objectives so
        this function will deal with that

        Parameters
        ----------
        funcs : dictionary of function values

        Returns
        -------
        obj : float or array
            Processed objective(s).
            """
        fobj = []
        for objKey in self.objectives.keys():
            if objKey in funcs:
                try: 
                    f = numpy.asscalar(numpy.squeeze(funcs[objKey]))
                except ValueError:
                    raise Error("The objective return value, '%s' must be a "
                                "scalar!"% objKey)
                # Store objective for printing later
                self.objectives[objKey].value = f
                if scaled:
                    f *= self.objectives[objKey].scale
                fobj.append(f)
            else:
                raise Error("The key for the objective, '%s' was not found." %
                            objKey)

        # Finally squeeze back out so we get a scalar for a single objective
        return numpy.squeeze(fobj)

    def processConstraints(self, fcon_in, scaled=True, dtype='d', natural=False):
        """
        **This function should not need to be called by the user**

        Parameters
        ----------
        fcon_in : dict
            Dictionary of constraint values

        scaled : bool
            Flag specifying if the returned array should be scaled by
            the pyOpt scaling. The only type this is not true is
            when the automatic derivatives are used

        dtype : str
            String specifiying the data type to return. Normally this
            is 'd' for a float. The complex-step derivative
            computations will call this function with 'D' to ensure
            that the complex peturbations pass through correctly.

        natural : bool
            Flag to specify if the data should be returned in the
            natural ordering. This is only used when computing
            gradient automatically with FD/CS.
            """

        if self.dummyConstraint:
            return numpy.array([0])
        
        # We REQUIRE that fcon_in is a dict:
        fcon = numpy.zeros(self.nCon, dtype=dtype)
        for iCon in self.constraints:
            con = self.constraints[iCon]
            if iCon in fcon_in:
                    
                # Make sure it is at least 1dimension:
                c = numpy.atleast_1d(fcon_in[iCon])
                if dtype == 'd':
                    c = numpy.real(c)
                # Make sure it is the correct size:
                if len(c) == self.constraints[iCon].ncon:
                    fcon[con.rs:con.re] = c
                else:
                    raise Error("%d constraint values were returned in "
                                "%s, but expected %d." % (
                                    len(fcon_in[iCon]), iCon,
                                    self.constraints[iCon].ncon))

                # Store constraint values for printing later
                con.value = copy.copy(c)
            else:
                raise Error("No constraint values were found for the "
                            "constraint '%s'."% iCon)

        # Perform scaling on the original jacobian:
        if scaled:
            fcon = self.conScale*fcon

        if natural:
            return fcon
        else:
            if self.nCon > 0:
                fcon = fcon[self.jacIndices]
                fcon = self.fact*fcon - self.offset
                return fcon
            else:
                return fcon
   
    def deProcessConstraints(self, fcon_in, scaled=True, dtype='d', natural=False):
        """
        **This function should not need to be called by the user**

        Parameters
        ----------
        fcon_in : array
            Array of constraint values to be converted into a dictionary

        scaled : bool
            Flag specifying if the returned array should be scaled by
            the pyOpt scaling. The only type this is not true is
            when the automatic derivatives are used

        dtype : str
            String specifiying the data type to return. Normally this
            is 'd' for a float. The complex-step derivative
            computations will call this function with 'D' to ensure
            that the complex peturbations pass through correctly.

        natural : bool
            Flag to specify if the input data is in the
            natural ordering. This is only used when computing
            gradient automatically with FD/CS.
            """

        if self.dummyConstraint:
            return {'dummy':0}

        if not natural:
            if self.nCon > 0:
                fcon_in += self.offset
                # Since self.fact elements are unit magnitude and the
                # values are either 1 or -1...
                fcon_in = self.fact * fcon_in
                # Undo the ordering
                fcon_in = fcon_in[self.jacIndicesInv]

        # Perform constraint scaling
        if scaled:
            fcon_in = fcon_in*self.conScale

        # We REQUIRE that fcon_in is an array:
        fcon = {}
        index = 0
        for iCon in self.constraints:
            con = self.constraints[iCon]
            fcon[iCon] = fcon_in[con.rs:con.re]

        return fcon

    def evaluateLinearConstraints(self, x, fcon):
        """
        This function is required for optimizers that do not explictly
        treat the linear constraints. For those optimizers, we will
        evaluate the linear constraints here. We place the values of
        the linear constraints in the fcon dictionary such that it
        appears as if the user evaluated these constraints. 

        Parameters
        ----------
        x : array
            This must be the processed x-vector from the optimizer

        fcon : dict
            Dictionary of the constraints. The linear constraints are
            to be added to this dictionary. 
            """

        # This is actually pretty easy; it's just a matvec with the
        # proper linearJacobian entry we've already computed
        for iCon in self.constraints:
            if self.constraints[iCon].linear:
                fcon[iCon] = multCOO(self.constraints[iCon].linearJacobian, x)

    def processObjectiveGradient(self, funcsSens):
        """
        **This function should not need to be called by the user**
        
        This generic function is used to assemble the objective
        gradient(s)

        Parameters
        ----------
        funcsSens : dict
            Dictionary of all function gradients. Just extract the
            objective(s) we need here. 
        """

        dvGroups = set(self.variables.keys())

        nobj = len(self.objectives)
        gobj = numpy.zeros((nobj, self.ndvs))

        iObj = 0
        for objKey in self.objectives.keys():
            if objKey in funcsSens:
                for dvGroup in funcsSens[objKey]:
                    if dvGroup in dvGroups:
                        # Now check that the array is the correct length:
                        ss = self.dvOffset[dvGroup]
                        tmp = numpy.array(funcsSens[objKey][dvGroup]).squeeze()
                        if tmp.size == ss[1]-ss[0]:
                            # Everything checks out so set:
                            gobj[iObj, ss[0]:ss[1]] = tmp * self.objectives[objKey].scale
                        else:
                            raise Error("The shape of the objective derivative "
                                        "for dvGroup '%s' is the incorrect "
                                        "length. Expecting a shape of %s but "
                                        "received a shape of %s."% (
                                            dvGroup, (ss[1]-ss[0],),
                                            funcsSens[objKey][dvGroup].shape))
                    else:
                        raise Error("The dvGroup key '%s' is not valid"% dvGroup)
            else:
                raise Error("The key for the objective gradient, '%s', was not found." %
                            objKey)
            iObj += 1

            # Note that we looped over the keys in funcsSens[objKey]
            # and not the variable keys since a variable key not in
            # funcsSens[objKey] will just be left to zero. We have
            # implictly assumed that the objective gradient is dense
            # and any keys that are provided are simply zero.
        # end (objective keys)

        # Do column scaling (dv scaling)
        gobj = self.invXScale * gobj

        # Finally squeeze back out so we get a 1D vector for a single objective
        return numpy.squeeze(gobj)

    def processConstraintJacobian(self, gcon):
        """
        **This function should not need to be called by the user**
        
        This generic function is used to assemble the entire
        constraint jacobian. The order of the constraint jacobian is
        in 'natural' ordering, that is the order the constraints have
        been added (mostly; since it can be different when constraints
        are added on different processors). 

        The input is gcon, which is dict or an array. The array format
        should only be used when the pyOpt_gradient class is used
        since this results in a desnse (and correctly oriented)
        jacobian. The user should NEVER return a desnse jacobian since
        this extremely fickle and easy to break. The dict 'gcon' must
        contain only the non-linear constraints jacobians; the linear
        ones will be added automatically.

        Parameters
        ----------
        gcon : array or dict
            Constraint gradients. Either a complete 2D array or a nested
            dictionary of gradients given with respect to the variables.

        Returns
        -------
        gcon : dict with csr data
            Return the jacobain in a sparse csr format. 
            can be easily converted to csc, coo or dense format as
            required by individual optimizers
            """

        # We don't have constraints at all! However we *may* have to
        # include a dummy constraint:
        if self.nCon == 0:
            if self.dummyConstraint:
                return convertToCSR(numpy.zeros((1, self.ndvs)))
            else:
                return numpy.zeros((0, self.ndvs), 'd')

        # For simplicity we just add the linear constraints into gcon
        # so they can be processed along with the rest:
        for iCon in self.constraints:
            if self.constraints[iCon].linear:
                gcon[iCon] = copy.deepcopy(self.constraints[iCon].jac)
                
        # We now know we must process as a dictionary. Below are the
        # lists for the matrix entris. 
        data = []
        row = []
        col = []
        ii = 0

        # Otherwise, process constraints in the dictionary form. 
        # Loop over all constraints:
        for iCon in self.constraints:
            con = self.constraints[iCon]

            if not con.name in gcon:
                raise Error("The jacobian for the constraint '%s' was "
                            "not found in the returned dictionary."%
                            con.name)

            if not con.partialReturnOk:
                # The keys in gcon[iCon] MUST match PRECISELY
                # the keys in con.wrt....The user told us they
                # would supply derivatives wrt to these sets, and
                # then didn't, so scold them. 
                for dvGroup in con.jac:
                    if dvGroup not in gcon[iCon]:
                        raise Error(
                            "Constraint '%s' was expecting a jacobain with "
                            "respect to dvGroup '%s' as was supplied in "
                            "addConGroup(). This was not found in the "
                            "constraint jacobian dictionary"% (
                                        con.name, dvGroup))
                    
            # Now loop over all required keys for this constraint:
            for dvGroup in con.wrt:
                # ss means 'start - stop'
                ss = self.dvOffset[dvGroup]
                ndvs = ss[1]-ss[0]
                if dvGroup in gcon[iCon]:
                    tmp = convertToCOO(gcon[iCon][dvGroup])
                else:
                    # This key is not returned. Just use the
                    # stored jacobian that contains zeros
                    tmp = con.jac[dvGroup]
                # Now check that the jacobian is the correct shape
                if not(tmp['shape'][0] == con.ncon and tmp['shape'][1] == ndvs):
                    raise Error("The shape of the supplied constraint "
                                "jacobian for constraint %s with respect to %s "
                                "is incorrect. "
                                "Expected an array of shape (%d, %d), but "
                                "received an array of shape (%d, %d)."% (
                                    con.name, dvGroup, con.ncon, ndvs, 
                                    tmp['shape'][0], tmp['shape'][1]))

                # Now check that supplied coo matrix has same length
                # of data array
                if len(tmp['coo'][2]) != len(con.jac[dvGroup]['coo'][2]):
                    raise Error("The number of nonzero elements for "
                                "constraint group '%s' with respect to %s "
                                "was not the correct size. The supplied "
                                "jacobian has %d nonzero "
                                "entries, but must contain %d nonzero "
                                "entries." % (con.name, dvGroup, len(tmp['coo'][2]),
                                              len(con.jac[dvGroup]['coo'][2])))

                # Include data from this jacobian chunk
                data.extend(tmp['coo'][IDATA])
                row.extend(tmp['coo'][IROW] + ii)
                col.extend(tmp['coo'][ICOL] + ss[0])
            # end for (dvGroup in constraint)
            ii += con.ncon
        # end for (constraint loop)

        # Finally, construct coo matrix, convert to csr and perform
        # row and column scaling.
        gcon = {'coo':[row, col, numpy.array(data)],
                'shape':[self.nCon, self.ndvs]}
        gcon = convertToCSR(gcon)
        scaleRows(gcon, self.conScale)
        scaleColumns(gcon, self.invXScale)

        return gcon

    def __str__(self):
        """
        Print Structured Optimization Problem
        """
        
        text = """\nOptimization Problem -- %s\n%s\n
        Objective Function: %s\n\n    Objectives:
        Name        Value        Optimum\n""" % (
        self.name, '='*80, self.objFun.__name__)

        for obj in self.objectives:
            lines = str(self.objectives[obj]).split('\n')
            text += lines[1] + '\n'

        text += """\n	Variables (c - continuous, i - integer, d - discrete):
           Name      Type       Value       Lower Bound  Upper Bound\n"""

        for dvGroup in self.variables:
            for var in self.variables[dvGroup]:
                lines = str(var).split('\n')
                text += lines[1] + '\n'

        print('	    Name        Type'+' '*25+'Bound\n'+'	 ')
        if len(self.constraints) > 0:
            text += """\n	Constraints (i - inequality, e - equality):
        Name    Type                    Bounds\n"""
            for iCon in self.constraints:
                text += str(self.constraints[iCon])

        return text
 
#==============================================================================
# Optimization Test
#==============================================================================
if __name__ == '__main__':
    
    print('Testing Optimization...')
    optprob = Optimization('Optimization Problem', {})
    
    
