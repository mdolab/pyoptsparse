# Standard Python modules
from collections import OrderedDict
import copy
import os
from typing import Callable, Dict, Iterable, List, Optional, Tuple, Union

# External modules
import numpy as np
from numpy import ndarray
from scipy.sparse import coo_matrix
from sqlitedict import SqliteDict

# Local modules
from .pyOpt_MPI import MPI
from .pyOpt_constraint import Constraint
from .pyOpt_objective import Objective
from .pyOpt_types import Dict1DType, Dict2DType, NumpyType
from .pyOpt_utils import (
    ICOL,
    IDATA,
    INFINITY,
    IROW,
    _broadcast_to_array,
    convertToCOO,
    convertToCSR,
    mapToCSR,
    scaleColumns,
    scaleRows,
)
from .pyOpt_variable import Variable


class Optimization:
    def __init__(self, name: str, objFun: Callable, comm=None, sens: Optional[Union[str, Callable]] = None):
        """
        The main purpose of this class is to describe the structure and
        potentially, sparsity pattern of an optimization problem.

        Parameters
        ----------
        name : str
            Name given to optimization problem.

        objFun : Python function handle
            Function handle used to evaluate the objective function.

        comm : MPI intra communication
            The communicator this problem will be solved on. This is
            required for both analysis when the objective is computed in
            parallel as well as to use the internal parallel gradient
            computations. Defaults to MPI.COMM_WORLD if not given.

        sens : str or python Function.
            Specify method to compute sensitivities.
        """
        self.name = name
        self.objFun = objFun
        self.sens = sens
        if comm is None:
            self.comm = MPI.COMM_WORLD
        else:
            self.comm = comm

        # Ordered dictionaries to keep track of variables and constraints
        self.variables: OrderedDict = OrderedDict()
        self.constraints: OrderedDict = OrderedDict()
        self.objectives: OrderedDict = OrderedDict()
        self.dvOffset: OrderedDict = OrderedDict()

        # Variables to be set in finalizeConstraints
        # have finalized the specification of the variable and the
        # constraints
        self.ndvs: int = 0
        self.conScale: ndarray = None
        self.nCon: int = 0
        self.nObj: int = 0
        self.invXScale: ndarray = None
        self.xOffset: ndarray = None
        self.dummyConstraint = False
        self.objectiveIdx: Dict[str, int] = {}
        self.finalized: bool = False
        self.jacIndices: ndarray = None
        self.fact: ndarray = None
        self.offset: ndarray = None

        # Store the Jacobian conversion maps
        self._jac_map_coo_to_csr = None

    def addVar(self, name: str, *args, **kwargs):
        """
        This is a convenience function. It simply calls addVarGroup()
        with nVars=1. Variables added with addVar() are returned as
        *scalars*.
        """
        self.addVarGroup(name, 1, *args, scalar=True, **kwargs)

    def checkVarName(self, varName: str) -> str:
        """
        Check if the desired variable name varName if has already been
        added. If it is has already been added, return a mangled name
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
            validName = f"{varName}_{i}"
            while validName in self.variables:
                i += 1
                validName = f"{varName}_{i}"
            return validName

    def checkConName(self, conName: str) -> str:
        """
        Check if the desired constraint name has already been
        added. If it is has already been added, return a mangled name
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
            validName = f"{conName}_{i}"
            while validName in self.constraints:
                i += 1
                validName = f"{conName}_{i}"
            return validName

    def addVarGroup(
        self,
        name: str,
        nVars: int,
        varType: str = "c",
        value=0.0,
        lower=None,
        upper=None,
        scale=1.0,
        offset=0.0,
        choices: List[str] = [],
        **kwargs,
    ):
        """
        Add a group of variables into a variable set. This is the main
        function used for adding variables to pyOptSparse.

        Parameters
        ----------
        name : str
            Name of variable group. This name should be unique across all the design variable groups

        nVars : int
            Number of design variables in this group.

        varType : str.
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
            necessary when design variables of widely varying
            magnitudes are used within the same
            optimization. Scalar/array usage is the same as value
            keyword.

        offset : scalar or array.  Define a user supplied offset
            variable for the design variable group.  This is often
            necessary when design variable has a large magnitude, but
            only changes a little about this value.

        choices : list
            Specify a list of choices for discrete design variables

        Examples
        --------
        >>> # Add a single design variable 'alpha'
        >>> optProb.addVar('alpha', varType='c', value=2.0, lower=0.0, upper=10.0, scale=0.1)
        >>> # Add 10 unscaled variables of 0.5 between 0 and 1 with name 'y'
        >>> optProb.addVarGroup('y', varType='c', value=0.5, lower=0.0, upper=1.0, scale=1.0)

        Notes
        -----
        Calling addVar() and addVarGroup(..., nVars=1, ...) are
        **NOT** equivalent! The variable added with addVar() will be
        returned as scalar, while variable returned from addVarGroup
        will be an array of length 1.

        It is recommended that the addVar() and addVarGroup() calls
        follow the examples above by including all the keyword
        arguments. This make it very clear the intent of the script's
        author. The type, value, lower, upper and scale should be
        given for all variables even if the default value is used.
        """
        self.finalized = False
        # Check that the nVars is > 0.
        if nVars < 1:
            raise ValueError(
                f"The 'nVars' argument to addVarGroup must be greater than or equal to 1. The bad DV is {name}."
            )

        # Check that the type is ok
        if varType not in ["c", "i", "d"]:
            raise ValueError("Type must be one of 'c' for continuous, 'i' for integer or 'd' for discrete.")

        value = _broadcast_to_array("value", value, nVars)
        lower = _broadcast_to_array("lower", lower, nVars, allow_none=True)
        upper = _broadcast_to_array("upper", upper, nVars, allow_none=True)
        scale = _broadcast_to_array("scale", scale, nVars)
        offset = _broadcast_to_array("offset", offset, nVars)

        # Determine if scalar i.e. it was called from addVar():
        scalar = kwargs.pop("scalar", False)

        # Now create all the variable objects
        varList = []
        for iVar in range(nVars):
            varName = f"{name}_{iVar}"
            varList.append(
                Variable(
                    varName,
                    varType=varType,
                    value=value[iVar],
                    lower=lower[iVar],
                    upper=upper[iVar],
                    scale=scale[iVar],
                    offset=offset[iVar],
                    scalar=scalar,
                    choices=choices,
                )
            )

        if name in self.variables:
            # Check that the variables happen to be the same
            if not len(self.variables[name]) == len(varList):
                raise KeyError(f"The supplied name '{name}' for a variable group has already been used!")
            for i in range(len(varList)):
                if not varList[i] == self.variables[name][i]:
                    raise KeyError(f"The supplied name '{name}' for a variable group has already been used!")
            # We we got here, we know that the variables we wanted to
            # add are **EXACTLY** the same so that's cool. We'll just
            # overwrite with the varList below.
        else:
            # Finally we set the variable list
            self.variables[name] = varList

    def delVar(self, name: str):
        """
        Delete a variable or variable group

        Parameters
        ----------
        name : str
           Name of variable or variable group to remove
        """
        self.finalized = False
        try:
            self.variables.pop(name)
        except KeyError:
            print(f"{name} was not a valid design variable name.")

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
                    if key not in procKeys:
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

    def addObj(self, name: str, *args, **kwargs):
        """
        Add Objective into Objectives Set
        """
        self.finalized = False
        self.objectives[name] = Objective(name, *args, **kwargs)

    def addCon(self, name: str, *args, **kwargs):
        """
        Convenience function. See addConGroup() for more information
        """
        self.addConGroup(name, 1, *args, **kwargs)

    def addConGroup(
        self,
        name: str,
        nCon: int,
        lower=None,
        upper=None,
        scale=1.0,
        linear: bool = False,
        wrt: Optional[Union[str, Iterable[str]]] = None,
        jac=None,
    ):
        r"""Add a group of variables into a variable set. This is the main
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
            advisable to have most optimization constraint around the
            same order of magnitude.

        linear : bool
            Flag to specify if this constraint is linear. If the
            constraint is linear, both the ``wrt`` and ``jac`` keyword
            arguments must be given to specify the constant portion of
            the constraint Jacobian.

            The intercept term of linear constraints must be supplied as
            part of the bound information. The linear constraint :math:`g_L \leq Ax + b \leq g_U`
            is equivalent to :math:`g_L - b \leq Ax \leq g_U - b`, and pyOptSparse requires
            the latter form. In this case, the arguments should be:

            .. code-block::

                jac = {"dvName" : A, ...}, lower = gL - b, upper = gU - b

        wrt : iterable (list, set, OrderedDict, array etc)
            'wrt' stand for stands for 'With Respect To'. This
            specifies for what dvs have non-zero Jacobian values
            for this set of constraints. The order is not important.

        jac : dictionary
            For linear and sparse non-linear constraints, the constraint
            Jacobian must be passed in. The structure of jac dictionary
            is as follows:

            .. code-block::

                {'dvName1':matrix1, 'dvName2', matrix1, ...}

            They keys of the Jacobian must correspond to the dvGroups
            given in the wrt keyword argument. The dimensions of each
            "chunk" of the constraint Jacobian must be consistent. For
            example, ``matrix1`` must have a shape of (nCon, nDvs) where
            nDVs is the total number of design variables in
            dvName1. ``matrix1`` may be a dense numpy array or it may be
            scipy sparse matrix. However, it is *HIGHLY* recommended
            that sparse constraints are supplied to pyOptSparse using
            the pyOptSparse's simplified sparse matrix format. The
            reason for this is that it is *impossible* for force scipy
            sparse matrices to keep a fixed sparsity pattern; if the
            sparsity pattern changes during an optimization, *IT WILL
            FAIL*.

            The three simplified pyOptSparse sparse matrix formats are
            summarized below:

            .. code-block::

                mat = {'coo':[row, col, data], 'shape':[nrow, ncols]} # A coo matrix
                mat = {'csr':[rowp, colind, data], 'shape':[nrow, ncols]} # A csr matrix
                mat = {'coo':[colp, rowind, data], 'shape':[nrow, ncols]} # A csc matrix

            Note that for nonlinear constraints (linear=False), the
            values themselves in the matrices in jac do not matter,
            but the sparsity structure **does** matter. It is
            imperative that entries that will at some point have
            non-zero entries have non-zero entries in jac
            argument. That is, we do not let the sparsity structure of
            the Jacobian change throughout the optimization. This
            stipulation is automatically checked internally.
        """
        self.finalized = False
        if name in self.constraints:
            raise KeyError(f"The supplied name '{name}' for a constraint group has already been used.")

        # Simply add constraint object
        self.constraints[name] = Constraint(name, nCon, linear, wrt, jac, lower, upper, scale)

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
        self.finalize()

        outDVs = {}
        for dvGroup in self.variables:
            nvar = len(self.variables[dvGroup])
            # If it is a single DV, return a scalar rather than a numpy array
            if nvar == 1:
                var = self.variables[dvGroup][0]
                outDVs[dvGroup] = var.value
            else:
                outDVs[dvGroup] = np.zeros(nvar)
                for i in range(nvar):
                    var = self.variables[dvGroup][i]
                    outDVs[dvGroup][i] = var.value
        # we convert the dict to array to scale everything consistently
        scaled_DV = self._mapXtoUser_Dict(outDVs)
        return scaled_DV

    def setDVs(self, inDVs):
        """
        Set one or more groups of design variables from a dictionary.
        In most common usage, this function is not required.

        Parameters
        ----------
        inDVs : dict
            The dictionary of variables. The keys are the names of the
            variable groups, and the values are the desired design
            variable values for each variable group.
        """
        self.finalize()
        x0 = self.getDVs()
        # overwrite subset of DVs with new values
        for dvGroup in inDVs:
            x0[dvGroup] = inDVs[dvGroup]
        # we process dicts to arrays to perform scaling in a uniform way
        # then process back to dict
        scaled_DV = self._mapXtoOpt_Dict(x0)
        for dvGroup in self.variables:
            if dvGroup in inDVs:
                nvar = len(self.variables[dvGroup])
                scalar = self.dvOffset[dvGroup][2]
                for i in range(nvar):
                    var = self.variables[dvGroup][i]
                    if scalar:
                        var.value = scaled_DV[dvGroup]
                    else:
                        # Must be an array
                        var.value = scaled_DV[dvGroup][i]

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
            hist = SqliteDict(histFile)
            if key is None:
                key = hist["last"]

            self.setDVs(hist[key]["xuser"])
            hist.close()
        else:
            raise FileNotFoundError(f"History file '{histFile}' not found!.")

    def printSparsity(self, verticalPrint=False):
        """
        This function prints an (ASCII) visualization of the Jacobian
        sparsity structure. This helps the user visualize what
        pyOptSparse has been given and helps ensure it is what the
        user expected. It is highly recommended this function be
        called before the start of every optimization to verify the
        optimization problem setup.

        Parameters
        ----------
        verticalPrint : bool
            True if the design variable names in the header should be printed
            vertically instead of horizontally. If true, this will make the
            constraint Jacobian print out more narrow and taller.

        Warnings
        --------
        This function is **collective** on the optProb comm. It is
        therefore necessary to call this function on **all**
        processors of the optProb comm.
        """
        self.finalize()

        if self.comm.rank != 0:
            return

        # Header describing what we are printing:
        print("+" + "-" * 78 + "-" + "+")
        print("|" + " " * 19 + "Sparsity structure of constraint Jacobian" + " " * 19 + "|")
        print("+" + "-" * 78 + "-" + "+")

        # We will do this with a 2d numpy array of characters since it
        # will make slicing easier

        # First determine the requried number of rows
        nRow = 1  # Header
        nRow += 1  # Line
        maxConNameLen = 0
        for iCon in self.constraints:
            nRow += 1  # Name
            con = self.constraints[iCon]
            maxConNameLen = max(maxConNameLen, len(con.name) + 6 + int(np.log10(con.ncon)) + 1)
            nRow += 1  # Line

        # And now the columns:
        nCol = maxConNameLen
        nCol += 2  # Space plus line
        varCenters = []
        longestNameLength = 0
        for dvGroup in self.variables:
            nvar = self.dvOffset[dvGroup][1] - self.dvOffset[dvGroup][0]

            # If printing vertically, put in a blank string of length 3
            if verticalPrint:
                var_str = "   "

            # Otherwise, put in the variable and its size
            else:
                var_str = f"{dvGroup} ({nvar})"

            # Find the length of the longest name for design variables
            longestNameLength = max(len(dvGroup), longestNameLength)

            varCenters.append(nCol + len(var_str) / 2 + 1)
            nCol += len(var_str)
            nCol += 2  # Spaces on either side
            nCol += 1  # Line

        txt = np.zeros((nRow, nCol), dtype=str)
        txt[:, :] = " "
        # Outline of the matrix on left and top
        txt[1, maxConNameLen + 1 : -1] = "-"
        txt[2:-1, maxConNameLen + 1] = "|"

        # Print the variable names:
        iCol = maxConNameLen + 2
        for dvGroup in self.variables:
            nvar = self.dvOffset[dvGroup][1] - self.dvOffset[dvGroup][0]
            if verticalPrint:
                var_str = "   "
            else:
                var_str = f"{dvGroup} ({nvar})"
            var_str_length = len(var_str)
            txt[0, iCol + 1 : iCol + var_str_length + 1] = list(var_str)
            txt[2:-1, iCol + var_str_length + 2] = "|"
            iCol += var_str_length + 3

        # Print the constraint names;
        iRow = 2

        for iCon in self.constraints:
            con = self.constraints[iCon]
            name = con.name
            if con.linear:
                name = name + "(L)"

            name = f"{name} ({con.ncon})"
            var_str_length = len(name)
            # The name
            txt[iRow, maxConNameLen - var_str_length : maxConNameLen] = list(name)

            # Now we write a 'X' if there is something there:
            varKeys = list(self.variables.keys())
            for dvGroup in range(len(varKeys)):
                if varKeys[dvGroup] in con.wrt:
                    txt[int(iRow), int(varCenters[dvGroup])] = "X"

            # The separator
            txt[iRow + 1, maxConNameLen + 1 :] = "-"
            iRow += 2

        # Corners - just to make it nice :-)
        txt[1, maxConNameLen + 1] = "+"
        txt[-1, maxConNameLen + 1] = "+"
        txt[1, -1] = "+"
        txt[-1, -1] = "+"

        # If we're printing vertically, add an additional text array on top
        # of the already created txt array
        if verticalPrint:
            # It has the same width and a height corresponding to the length
            # of the longest design variable name
            newTxt = np.zeros((longestNameLength + 1, nCol), dtype=str)
            newTxt[:, :] = " "
            txt = np.vstack((newTxt, txt))

            # Loop through the letters in the longest design variable name
            # and add the letters for each design variable
            for i in range(longestNameLength + 2):
                # Make a space between the name and the size
                if i >= longestNameLength:
                    txt[i, :] = " "

                # Loop through each design variable
                for j, dvGroup in enumerate(self.variables):
                    # Print a letter in the name if any remain
                    if i < longestNameLength and i < len(dvGroup):
                        txt[i, int(varCenters[j])] = dvGroup[i]

                    # Format and print the size of the design variable
                    elif i > longestNameLength:
                        var_str = "(" + str(self.dvOffset[dvGroup][1] - self.dvOffset[dvGroup][0]) + ")"
                        half_length = len(var_str) / 2
                        k = int(varCenters[j])
                        txt[i, int(k - half_length + 1) : int(k - half_length + 1 + len(var_str))] = list(var_str)

        for i in range(len(txt)):
            print("".join(txt[i]))

    def getDVConIndex(self, startIndex: int = 1, printIndex: bool = True) -> Tuple[OrderedDict, OrderedDict]:
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
                dvIndex[dvGroup] = [ind0, ind1 - 1]

        # Get the begin and end index (inclusive) of constraints
        conIndex = OrderedDict()
        conCounter = startIndex
        for iCon in self.constraints:
            n = self.constraints[iCon].ncon
            if n == 1:
                conIndex[iCon] = [conCounter]
            else:
                conIndex[iCon] = [conCounter, conCounter + n - 1]
            conCounter += n

        # Print them all to terminal
        if printIndex and self.comm.rank == 0:
            print("### DESIGN VARIABLES ###")
            for dvGroup in dvIndex:
                print(dvGroup, dvIndex[dvGroup])
            print("### CONSTRAINTS ###")
            for conKey in conIndex:
                print(conKey, conIndex[conKey])

        return dvIndex, conIndex

    # =======================================================================
    #       All the functions from here down should not need to be called
    #       by the user. Most functions are public since the individual
    #       optimizers need to be able to call them
    # =======================================================================

    def finalize(self):
        """
        This is a helper function which will only finalize the optProb if it's not already finalized.
        """
        if not self.finalized:
            self._finalizeDesignVariables()
            self._finalizeConstraints()
            self.finalized = True

    def _finalizeDesignVariables(self):
        """
        Communicate design variables potentially from different
        processors and form the DVOffset dict.

        Warnings
        --------
        This should not be called directly. Instead, call self.finalize()
        to ensure that both design variables and constraints are properly finalized.
        """

        # First thing we need is to determine the consistent set of
        # variables from all processors.
        self.variables = self._reduceDict(self.variables)

        dvCounter = 0
        self.dvOffset = OrderedDict()

        for dvGroup in self.variables:
            n = len(self.variables[dvGroup])
            self.dvOffset[dvGroup] = [dvCounter, dvCounter + n, self.variables[dvGroup][0].scalar]
            dvCounter += n
        self.ndvs = dvCounter

    def _finalizeConstraints(self):
        """
        There are several functions for this routine:

        1. Determine the number of constraints
        2. Determine the final scaling array for the design variables
        3. Determine if it is possible to return a complete dense
           Jacobian. Most of this time, we should be using the dictionary-
           based return

        Warnings
        --------
        This should not be called directly. Instead, call self.finalize()
        to ensure that both design variables and constraints are properly finalized.
        """
        # reset these counters
        self.nObj = 0
        self.nCon = 0
        # First thing we need is to determine the consistent set of
        # constraints from all processors
        self.constraints = self._reduceDict(self.constraints)

        # ----------------------------------------------------
        # Step 1. Determine number of constraints and scaling:
        # ----------------------------------------------------

        # Determine number of constraints
        for iCon in self.constraints:
            self.nCon += self.constraints[iCon].ncon

        # Loop over the constraints assigning the row start (rs) and
        # row end (re) values. The actual ordering depends on if
        # constraints are reordered or not.
        rowCounter = 0
        conScale = np.zeros(self.nCon)
        for iCon in self.constraints:
            con = self.constraints[iCon]
            con.finalize(self.variables, self.dvOffset, rowCounter)
            rowCounter += con.ncon
            conScale[con.rs : con.re] = con.scale

        if self.nCon > 0:
            self.conScale = conScale
        else:
            self.conScale = None

        # -----------------------------------------
        # Step 2a. Assemble design variable scaling
        # -----------------------------------------
        xscale = []
        for dvGroup in self.variables:
            for var in self.variables[dvGroup]:
                xscale.append(var.scale)
        self.invXScale = 1.0 / np.array(xscale)

        # -----------------------------------------
        # Step 2a. Assemble design variable offset
        # -----------------------------------------
        xoffset = []
        for dvGroup in self.variables:
            for var in self.variables[dvGroup]:
                xoffset.append(var.offset)
        self.xOffset = np.array(xoffset)

        # --------------------------------------
        # Step 3. Map objective names to indices
        # --------------------------------------
        for idx, objKey in enumerate(self.objectives):
            self.objectiveIdx[objKey] = idx
            self.nObj += 1

        # ---------------------------------------------
        # Step 4. Final Jacobian for linear constraints
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

                    row.extend(con.jac[dvGroup]["coo"][IROW])
                    col.extend(con.jac[dvGroup]["coo"][ICOL] + ss[0])
                    data.extend(con.jac[dvGroup]["coo"][IDATA])

                # Now create a coo, convert to CSR and store
                con.linearJacobian = coo_matrix((data, (row, col)), shape=[con.ncon, self.ndvs]).tocsr()

    def getOrdering(
        self, conOrder: List[str], oneSided: bool, noEquality: bool = False
    ) -> Tuple[ndarray, ndarray, ndarray, ndarray]:
        """
        Internal function that is used to produce a index list that
        reorders the constraints the way a particular optimizer needs.

        Parameters
        ----------
        conOrder : list
            This must contain the following 4 strings: 'ni', 'li',
            'ne', 'le' which stand for nonlinear inequality, linear
            inequality, nonlinear equality and linear equality. This
            defines the order that the optimizer wants the constraints

        oneSided : bool
           Flag to do all constraints as one-sided instead of two
           sided. Most optimizers need this but some can deal with the
           two-sided constraints properly (snopt and ipopt for
           example)

        noEquality : bool
           Flag to split equality constraints into two inequality
           constraints. Some optimizers (CONMIN for example) can't do
           equality constraints explicitly.
        """

        # Now for the fun part determine what *actual* order the
        # constraints need to be in: We recognize the following
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
                return np.array([], "d")

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

                if conType == "ne" and not con.linear:
                    if noEquality:
                        # Expand Equality constraint to two:
                        indices.extend(con.rs + econ["ind"])
                        fact.extend(econ["fact"])
                        lower.extend(econ["value"])
                        upper.extend(econ["value"])
                        # ....And the other side
                        indices.extend(con.rs + econ["ind"])
                        fact.extend(-1.0 * econ["fact"])
                        lower.extend(econ["value"])
                        upper.extend(econ["value"])

                    else:
                        indices.extend(con.rs + econ["ind"])
                        fact.extend(econ["fact"])
                        lower.extend(econ["value"])
                        upper.extend(econ["value"])

                if conType == "ni" and not con.linear:
                    indices.extend(con.rs + icon["ind"])
                    fact.extend(icon["fact"])
                    lower.extend(icon["lower"])
                    upper.extend(icon["upper"])

                if conType == "le" and con.linear:
                    if noEquality:
                        # Expand Equality constraint to two:
                        indices.extend(con.rs + econ["ind"])
                        fact.extend(econ["fact"])
                        lower.extend([-INFINITY] * len(econ["fact"]))
                        upper.extend(econ["value"])
                        # ....And the other side
                        indices.extend(con.rs + econ["ind"])
                        fact.extend(-1.0 * econ["fact"])
                        lower.extend([-INFINITY] * len(econ["fact"]))
                        upper.extend(-econ["value"])
                    else:
                        indices.extend(con.rs + econ["ind"])
                        fact.extend(econ["fact"])
                        lower.extend(econ["value"])
                        upper.extend(econ["value"])

                if conType == "li" and con.linear:
                    indices.extend(con.rs + icon["ind"])
                    fact.extend(icon["fact"])
                    lower.extend(icon["lower"])
                    upper.extend(icon["upper"])

        return np.array(indices), np.array(lower), np.array(upper), np.array(fact)

    def processXtoDict(self, x: ndarray) -> OrderedDict:
        """
        Take the flattened array of variables in 'x' and return a
        dictionary of variables keyed on the name of each variable.

        Parameters
        ----------
        x : array
            Flattened array from optimizer

        Warnings
        --------
        This function should not need to be called by the user
        """
        xg = OrderedDict()
        imax = 0
        for dvGroup in self.variables:
            istart = self.dvOffset[dvGroup][0]
            iend = self.dvOffset[dvGroup][1]
            scalar = self.dvOffset[dvGroup][2]
            imax = max(imax, iend)
            try:
                if scalar:
                    xg[dvGroup] = x[..., istart]
                else:
                    xg[dvGroup] = x[..., istart:iend].copy()
            except IndexError as e:
                raise ValueError("Error processing x. There is a mismatch in the number of variables.") from e
        if imax != self.ndvs:
            raise ValueError("Error processing x. There is a mismatch in the number of variables.")
        return xg

    def processXtoVec(self, x: dict) -> ndarray:
        """
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

        Warnings
        --------
        This function should not need to be called by the user
        """
        x_array = np.zeros(self.ndvs)
        imax = 0
        for dvGroup in self.variables:
            istart = self.dvOffset[dvGroup][0]
            iend = self.dvOffset[dvGroup][1]
            imax = max(imax, iend)
            scalar = self.dvOffset[dvGroup][2]
            try:
                if scalar:
                    x_array[..., istart] = x[dvGroup]
                else:
                    x_array[..., istart:iend] = x[dvGroup]
            except IndexError as e:
                raise ValueError("Error deprocessing x. There is a mismatch in the number of variables.") from e
        if imax != self.ndvs:
            raise ValueError("Error deprocessing x. There is a mismatch in the number of variables.")

        return x_array

    def processObjtoVec(self, funcs: Dict1DType, scaled: bool = True) -> NumpyType:
        """
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

        Warnings
        --------
        This function should not need to be called by the user
        """
        fobj = []
        for objKey in self.objectives.keys():
            if objKey in funcs:
                try:
                    f = np.squeeze(funcs[objKey]).item()
                except ValueError as e:
                    raise ValueError(f"The objective return value, '{objKey}' must be a scalar!") from e
                # Store objective for printing later
                self.objectives[objKey].value = np.real(f)
                fobj.append(f)
            else:
                raise KeyError(f"The key for the objective, '{objKey}' was not found.")

        # scale the objective
        if scaled:
            fobj = self._mapObjtoOpt(fobj)
        # Finally squeeze back out so we get a scalar for a single objective
        return np.squeeze(fobj)

    def processObjtoDict(self, fobj_in: NumpyType, scaled: bool = True) -> Dict1DType:
        """
        This function converts the objective in array form
        to the corresponding dictionary form.

        Parameters
        ----------
        fobj_in : float or ndarray
            The objective in array format. In the case of a single objective,
            a float can also be accepted.
        scaled : bool
            Flag specifying if the returned dictionary should be scaled by
            the pyOpt scaling.

        Returns
        -------
        fobj : dictionary
            The dictionary form of fobj_in, which is just a key:value pair
            for each objective.
        """
        fobj = {}
        fobj_in = np.atleast_1d(fobj_in)
        for objKey in self.objectives.keys():
            iObj = self.objectiveIdx[objKey]
            try:
                fobj[objKey] = fobj_in[iObj]
            except IndexError as e:
                raise ValueError("The input array shape is incorrect!") from e
        if scaled:
            fobj = self._mapObjtoOpt(fobj)
        return fobj

    def processContoVec(
        self, fcon_in: Dict1DType, scaled: bool = True, dtype: str = "d", natural: bool = False
    ) -> ndarray:
        """A function that converts a dictionary of constraints into a vector

        Parameters
        ----------
        fcon_in : dict
            Dictionary of constraint values

        scaled : bool
            Flag specifying if the returned array should be scaled by
            the pyOpt scaling. The only type this is not true is
            when the automatic derivatives are used

        dtype : str
            String specifying the data type to return. Normally this
            is 'd' for a float. The complex-step derivative
            computations will call this function with 'D' to ensure
            that the complex perturbations pass through correctly.

        natural : bool
            Flag to specify if the data should be returned in the
            natural ordering. This is only used when computing
            gradient automatically with FD/CS.

        Warnings
        --------
        This function should not need to be called by the user
        """

        if self.dummyConstraint:
            return np.array([0])

        # We REQUIRE that fcon_in is a dict:
        fcon = np.zeros(self.nCon, dtype=dtype)
        for iCon in self.constraints:
            con = self.constraints[iCon]
            if iCon in fcon_in:
                # Make sure it is at least 1-dimensional:
                c = np.atleast_1d(fcon_in[iCon])
                if dtype == "d":
                    c = np.real(c)
                # Make sure it is the correct size:
                if c.shape[-1] == self.constraints[iCon].ncon:
                    fcon[..., con.rs : con.re] = c
                else:
                    raise ValueError(
                        f"{len(fcon_in[iCon])} constraint values were returned in {iCon}, "
                        + f"but expected {self.constraints[iCon].ncon}."
                    )

                # Store constraint values for printing later
                con.value = np.real(copy.copy(c))
            else:
                raise KeyError(f"No constraint values were found for the constraint '{iCon}'.")

        # Perform scaling on the original Jacobian:
        if scaled:
            fcon = self._mapContoOpt(fcon)

        if natural:
            return fcon
        else:
            if self.nCon > 0:
                fcon = fcon[..., self.jacIndices]
                fcon = self.fact * fcon - self.offset
                return fcon
            else:
                return fcon

    def processContoDict(
        self, fcon_in: ndarray, scaled: bool = True, dtype: str = "d", natural: bool = False, multipliers: bool = False
    ) -> Dict1DType:
        """A function that converts an array of constraints into a dictionary

        Parameters
        ----------
        fcon_in : array
            Array of constraint values to be converted into a dictionary

        scaled : bool
            Flag specifying if the returned array should be scaled by
            the pyOpt scaling. The only time this is not true is
            when the automatic derivatives are used

        dtype : str
            String specifying the data type to return. Normally this
            is 'd' for a float. The complex-step derivative
            computations will call this function with 'D' to ensure
            that the complex perturbations pass through correctly.

        natural : bool
            Flag to specify if the input data is in the
            natural ordering. This is only used when computing
            gradient automatically with FD/CS.

        multipliers : bool
            Flag that indicates whether this deprocessing is for the
            multipliers or the constraint values. In the case of multipliers,
            no constraint offset should be applied.

        Warnings
        --------
        This function should not need to be called by the user
        """

        if self.dummyConstraint:
            return {"dummy": 0}

        if not hasattr(self, "jacIndicesInv"):
            self.jacIndicesInv = np.argsort(self.jacIndices)

        # Unscale the nonlinear constraints
        if not natural:
            if self.nCon > 0:
                m = len(self.jacIndices)
                # Apply the offset (if this is for constraint values)
                if not multipliers:
                    fcon_in[:m] += self.offset

                # Since self.fact elements are unit magnitude and the
                # values are either 1 or -1...
                fcon_in[:m] = self.fact * fcon_in[:m]

        # Perform constraint scaling
        if scaled:
            m = len(self.jacIndices)
            fcon_in[:m] = fcon_in[:m] * self.conScale[self.jacIndices]

        fcon_unique = fcon_in
        if multipliers:
            fcon_unique = np.zeros(self.nCon)
            for i, j in enumerate(self.jacIndices):
                if np.abs(fcon_unique[j]) < np.abs(fcon_in[i]):
                    fcon_unique[j] = fcon_in[i]

        # We REQUIRE that fcon_in is an array:
        fcon = {}
        for iCon in self.constraints:
            con = self.constraints[iCon]
            fcon[iCon] = fcon_unique[..., con.rs : con.re]

        return fcon

    def evaluateLinearConstraints(self, x: ndarray, fcon: Dict1DType):
        """
        This function is required for optimizers that do not explicitly
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
                fcon[iCon] = self.constraints[iCon].linearJacobian.dot(x)

    def processObjectiveGradient(self, funcsSens: Dict2DType) -> NumpyType:
        """
        This generic function is used to assemble the objective
        gradient(s)

        Parameters
        ----------
        funcsSens : dict
            Dictionary of all function gradients. Just extract the
            objective(s) we need here.

        Warnings
        --------
        This function should not need to be called by the user
        """

        dvGroups = set(self.variables.keys())
        gobj = np.zeros((self.nObj, self.ndvs))

        iObj = 0
        for objKey in self.objectives.keys():
            if objKey in funcsSens:
                for dvGroup in funcsSens[objKey]:
                    if dvGroup in dvGroups:
                        # Now check that the array is the correct length:
                        ss = self.dvOffset[dvGroup]
                        tmp = np.array(funcsSens[objKey][dvGroup]).squeeze()
                        if tmp.size == ss[1] - ss[0]:
                            # Everything checks out so set:
                            gobj[iObj, ss[0] : ss[1]] = tmp
                        else:
                            raise ValueError(
                                f"The shape of the objective derivative for dvGroup '{dvGroup}' is the incorrect length. "
                                + f"Expecting a shape of {(ss[1] - ss[0],)} but received a shape of {funcsSens[objKey][dvGroup].shape}."
                            )
                    else:
                        raise KeyError(f"The dvGroup key '{dvGroup}' is not valid")
            else:
                raise KeyError(f"The key for the objective gradient, '{objKey}', was not found.")
            iObj += 1

        # Note that we looped over the keys in funcsSens[objKey]
        # and not the variable keys since a variable key not in
        # funcsSens[objKey] will just be left to zero. We have
        # implicitly assumed that the objective gradient is dense
        # and any keys that are provided are simply zero.
        # end (objective keys)

        # Do scaling
        gobj = self._mapObjGradtoOpt(gobj)

        # Finally squeeze back out so we get a 1D vector for a single objective
        return np.squeeze(gobj)

    def processConstraintJacobian(self, gcon):
        """
        This generic function is used to assemble the entire
        constraint Jacobian. The order of the constraint Jacobian is
        in 'natural' ordering, that is the order the constraints have
        been added (mostly; since it can be different when constraints
        are added on different processors).

        The input is gcon, which is dict or an array. The array format
        should only be used when the pyOpt_gradient class is used
        since this results in a dense (and correctly oriented)
        Jacobian. The user should NEVER return a dense Jacobian since
        this extremely fickle and easy to break. The dict 'gcon' must
        contain only the non-linear constraints Jacobians; the linear
        ones will be added automatically.

        Parameters
        ----------
        gcon : array or dict
            Constraint gradients. Either a complete 2D array or a nested
            dictionary of gradients given with respect to the variables.

        Returns
        -------
        gcon : dict with csr data
            Return the Jacobian in a sparse csr format.
            can be easily converted to csc, coo or dense format as
            required by individual optimizers

        Warnings
        --------
        This function should not need to be called by the user
        """

        # We don't have constraints at all! However we *may* have to
        # include a dummy constraint:
        if self.nCon == 0:
            if self.dummyConstraint:
                return convertToCSR(np.zeros((1, self.ndvs)))
            else:
                return np.zeros((0, self.ndvs), "d")

        # For simplicity we just add the linear constraints into gcon
        # so they can be processed along with the rest:
        for iCon in self.constraints:
            if self.constraints[iCon].linear:
                gcon[iCon] = copy.deepcopy(self.constraints[iCon].jac)

        # We now know we must process as a dictionary. Below are the
        # lists for the matrix entries.
        data = []
        row = []
        col = []
        ii = 0

        # Otherwise, process constraints in the dictionary form.
        # Loop over all constraints:
        for iCon in self.constraints:
            con = self.constraints[iCon]

            # Now loop over all required keys for this constraint:
            for dvGroup in con.wrt:
                # ss means 'start - stop'
                ss = self.dvOffset[dvGroup]
                ndvs = ss[1] - ss[0]

                gotDerivative = False
                if dvGroup in gcon[iCon]:
                    tmp = convertToCOO(gcon[iCon][dvGroup])
                    gotDerivative = True
                else:
                    raise KeyError(
                        f"The constraint Jacobian entry for '{con.name}' with respect to '{dvGroup}', "
                        + "as was defined in addConGroup(), was not found in constraint Jacobian dictionary provided."
                    )
                if not gotDerivative:
                    # All keys for this constraint must be returned
                    # since the user has explicitly specified the wrt.
                    if not con.partialReturnOk:
                        raise ValueError(
                            f"Constraint '{con.name}' was expecting a Jacobian with respect to dvGroup "
                            + f"'{dvGroup}' as was supplied in addConGroup(). "
                            + "This was not found in the constraint Jacobian dictionary"
                        )
                    else:
                        # This key is not returned. Just use the
                        # stored Jacobian that contains zeros
                        tmp = con.jac[dvGroup]

                # Now check that the Jacobian is the correct shape
                if not (tmp["shape"][0] == con.ncon and tmp["shape"][1] == ndvs):
                    raise ValueError(
                        f"The shape of the supplied constraint Jacobian for constraint {con.name} with respect to {dvGroup} is incorrect. "
                        + f"Expected an array of shape ({con.ncon}, {ndvs}), but received an array of shape ({tmp['shape'][0]}, {tmp['shape'][1]})."
                    )

                # Now check that supplied coo matrix has same length
                # of data array
                if len(tmp["coo"][2]) != len(con.jac[dvGroup]["coo"][2]):
                    raise ValueError(
                        f"The number of nonzero elements for constraint group '{con.name}' with respect to {dvGroup} was not the correct size. "
                        + f"The supplied Jacobian has {len(tmp['coo'][2])} nonzero entries, but must contain {len(con.jac[dvGroup]['coo'][2])} nonzero entries."
                    )

                # Include data from this Jacobian chunk
                data.append(tmp["coo"][IDATA])
                row.append(tmp["coo"][IROW] + ii)
                col.append(tmp["coo"][ICOL] + ss[0])
            # end for (dvGroup in constraint)
            ii += con.ncon
        # end for (constraint loop)

        # now flatten all the data into a single array
        data = np.concatenate(data).ravel()
        row = np.concatenate(row).ravel()
        col = np.concatenate(col).ravel()

        # Finally, construct CSR matrix from COO data and perform
        # row and column scaling.
        if self._jac_map_coo_to_csr is None:
            gcon = {"coo": [row, col, np.array(data)], "shape": [self.nCon, self.ndvs]}
            self._jac_map_coo_to_csr = mapToCSR(gcon)

        gcon = {
            "csr": (
                self._jac_map_coo_to_csr[IROW],
                self._jac_map_coo_to_csr[ICOL],
                np.array(data)[self._jac_map_coo_to_csr[IDATA]],
            ),
            "shape": [self.nCon, self.ndvs],
        }

        self._mapConJactoOpt(gcon)

        return gcon

    def _mapObjGradtoOpt(self, gobj: ndarray) -> ndarray:
        gobj_return = np.copy(gobj)
        for objKey in self.objectives:
            iObj = self.objectiveIdx[objKey]
            gobj_return[iObj, :] *= self.objectives[objKey].scale
        gobj_return *= self.invXScale
        return gobj_return

    def _mapContoOpt(self, fcon: ndarray) -> ndarray:
        return fcon * self.conScale

    def _mapContoUser(self, fcon: ndarray) -> ndarray:
        return fcon / self.conScale

    def _mapObjtoOpt(self, fobj: ndarray) -> ndarray:
        fobj_return = np.copy(np.atleast_1d(fobj))
        for objKey in self.objectives:
            iObj = self.objectiveIdx[objKey]
            fobj_return[iObj] *= self.objectives[objKey].scale
        return fobj_return

    def _mapObjtoUser(self, fobj: ndarray) -> ndarray:
        fobj_return = np.copy(np.atleast_1d(fobj))
        for objKey in self.objectives:
            iObj = self.objectiveIdx[objKey]
            fobj_return[iObj] /= self.objectives[objKey].scale
        return fobj_return

    def _mapConJactoOpt(self, gcon: ndarray) -> ndarray:
        """
        The mapping is done in memory, without any return.
        """
        scaleRows(gcon, self.conScale)
        scaleColumns(gcon, self.invXScale)

    def _mapConJactoUser(self, gcon: ndarray) -> ndarray:
        """
        The mapping is done in memory, without any return.
        """
        scaleRows(gcon, 1 / self.conScale)
        scaleColumns(gcon, 1 / self.invXScale)

    def _mapXtoOpt(self, x: ndarray) -> ndarray:
        """
        This performs the user-space to optimizer mapping for the DVs.
        All inputs/outputs are numpy arrays.
        """
        return (x - self.xOffset) / self.invXScale

    def _mapXtoUser(self, x: ndarray) -> ndarray:
        """
        This performs the optimizer to user-space mapping for the DVs.
        All inputs/outputs are numpy arrays.
        """
        return x * self.invXScale + self.xOffset

    # these are the dictionary-based versions of the mapping functions
    def _mapXtoUser_Dict(self, xDict: Dict1DType) -> Dict1DType:
        x = self.processXtoVec(xDict)
        x_user = self._mapXtoUser(x)
        return self.processXtoDict(x_user)

    def _mapXtoOpt_Dict(self, xDict: Dict1DType) -> Dict1DType:
        x = self.processXtoVec(xDict)
        x_opt = self._mapXtoOpt(x)
        return self.processXtoDict(x_opt)

    def _mapObjtoUser_Dict(self, objDict: Dict1DType) -> Dict1DType:
        obj = self.processObjtoVec(objDict, scaled=False)
        obj_user = self._mapObjtoUser(obj)
        return self.processObjtoDict(obj_user, scaled=False)

    def _mapObjtoOpt_Dict(self, objDict: Dict1DType) -> Dict1DType:
        obj = self.processObjtoVec(objDict, scaled=False)
        obj_opt = self._mapObjtoOpt(obj)
        return self.processObjtoDict(obj_opt, scaled=False)

    def _mapContoUser_Dict(self, conDict: Dict1DType) -> Dict1DType:
        con = self.processContoVec(conDict, scaled=False, natural=True)
        con_user = self._mapContoUser(con)
        return self.processContoDict(con_user, scaled=False, natural=True)

    def _mapContoOpt_Dict(self, conDict: Dict1DType) -> Dict1DType:
        con = self.processContoVec(conDict, scaled=False, natural=True)
        con_opt = self._mapContoOpt(con)
        return self.processContoDict(con_opt, scaled=False, natural=True)

    def summary_str(self, minimal_print=False, print_multipliers=False):
        """
        Print Structured Optimization Problem

        Parameters
        ----------
        minimal_print : bool
            Flag to specify if the printed results should only include
            variables and constraints with a non-empty status
            (for example a violated bound).
            This defaults to False, which will print all results.
        print_multipliers : bool
            If True, print the Lagrange multipliers associated with the constraints.
        """
        TOL = 1.0e-6

        text = (
            f"\n\nOptimization Problem -- {self.name}\n{'=' * 80}\n    Objective Function: {self.objFun.__name__}\n\n"
        )
        text += "\n   Objectives\n"

        num_c = max(len(obj) for obj in self.objectives)
        fmt = "    {0:>7s}  {1:{width}s}   {2:>14s}\n"
        text += fmt.format("Index", "Name", "Value", width=num_c)
        fmt = "    {0:>7d}  {1:{width}s}   {2:>14.6E}\n"
        for idx, name in enumerate(self.objectives):
            obj = self.objectives[name]
            text += fmt.format(idx, obj.name, obj.value, width=num_c)

        # Find the longest name in the variables
        num_c = 0
        for varname in self.variables:
            for var in self.variables[varname]:
                num_c = max(len(var.name), num_c)

        fmt = "    {0:>7s}  {1:{width}s}   {2:>4s}   {3:>14}   {4:>14}   {5:>14}   {6:>8s}\n"
        text += "\n   Variables (c - continuous, i - integer, d - discrete)\n"
        text += fmt.format("Index", "Name", "Type", "Lower Bound", "Value", "Upper Bound", "Status", width=num_c)
        fmt = "    {0:7d}  {1:{width}s}   {2:>4s}   {3:14.6E}   {4:14.6E}   {5:14.6E}   {6:>8s}\n"

        idx = 0
        for varname in self.variables:
            for var in self.variables[varname]:
                if var.type in ["c", "i"]:
                    value = var.value
                    lower = var.lower if var.lower is not None else -1.0e20
                    upper = var.upper if var.upper is not None else 1.0e20
                    status = ""
                    dL = value - lower
                    if dL > TOL:
                        pass
                    elif dL < -TOL:
                        # In violation of lower bound
                        status += "L"
                    else:
                        # Active lower bound
                        status += "l"
                    dU = upper - value
                    if dU > TOL:
                        pass
                    elif dU < -TOL:
                        # In violation of upper bound
                        status += "U"
                    else:
                        # Active upper bound
                        status += "u"
                elif var.type == "d":
                    choices = var.choices
                    value = choices[int(var.value)]
                    lower = min(choices)
                    upper = max(choices)
                    status = ""
                else:
                    raise ValueError(f"Unrecognized type for variable {var.name}: {var.type}")

                if not minimal_print or status:
                    text += fmt.format(idx, var.name, var.type, lower, value, upper, status, width=num_c)
                idx += 1

        if len(self.constraints) > 0:
            # must be an instance of the Solution class
            if print_multipliers and self.lambdaStar is not None:
                lambdaStar = self.lambdaStar
                lambdaStar_label = "Lagrange Multiplier"
            else:
                # the optimizer did not set the lagrange multipliers so set them to something obviously wrong
                lambdaStar = {}
                for c in self.constraints:
                    lambdaStar[c] = [9e100] * self.constraints[c].ncon
                lambdaStar_label = "Lagrange Multiplier (N/A)"

            text += "\n   Constraints (i - inequality, e - equality)\n"
            # Find the longest name in the constraints
            num_c = max(len(self.constraints[i].name) for i in self.constraints)
            fmt = "    {0:>7s}  {1:{width}s} {2:>4s} {3:>14}  {4:>14}  {5:>14}  {6:>8s}  {7:>14s}\n"
            text += fmt.format(
                "Index", "Name", "Type", "Lower", "Value", "Upper", "Status", lambdaStar_label, width=num_c
            )
            fmt = "    {0:7d}  {1:{width}s} {2:>4s} {3:>14.6E}  {4:>14.6E}  {5:>14.6E}  {6:>8s}  {7:>14.5E}\n"
            idx = 0
            for con_name in self.constraints:
                c = self.constraints[con_name]
                for j in range(c.ncon):
                    lower = c.lower[j] if c.lower[j] is not None else -1.0e20
                    upper = c.upper[j] if c.upper[j] is not None else 1.0e20
                    value = c.value[j]
                    status = ""
                    typ = "e" if j in c.equalityConstraints["ind"] else "i"
                    if typ == "e":
                        if abs(value - upper) > TOL:
                            status = "E"
                    else:
                        dL = value - lower
                        if dL > TOL:
                            pass
                        elif dL < -TOL:
                            # In violation of lower bound
                            status += "L"
                        else:
                            # Active lower bound
                            status += "l"

                        dU = upper - value
                        if dU > TOL:
                            pass
                        elif dU < -TOL:
                            # In violation of upper bound
                            status += "U"
                        else:
                            # Active upper bound
                            status += "u"

                    if not minimal_print or status:
                        text += fmt.format(
                            idx, c.name, typ, lower, value, upper, status, lambdaStar[con_name][j], width=num_c
                        )
                    idx += 1

        return text

    def __str__(self):
        return self.summary_str(minimal_print=False, print_multipliers=False)

    def __getstate__(self) -> dict:
        """
        This is used for serializing class instances.
        The un-serializable fields are deleted first.
        """
        d = copy.copy(self.__dict__)
        for key in ["comm"]:
            if key in d.keys():
                del d[key]
        return d
