# Standard Python modules
from collections import OrderedDict
import copy
import datetime
from enum import Enum
import os
import shutil
import tempfile
import time
from typing import Any, Callable, Dict, List, Optional, Union

# External modules
from baseclasses import BaseSolver
import numpy as np
from numpy import ndarray

# Local modules
from .pyOpt_MPI import MPI
from .pyOpt_error import pyOptSparseWarning
from .pyOpt_gradient import Gradient
from .pyOpt_history import History
from .pyOpt_optimization import Optimization
from .pyOpt_solution import Solution
from .pyOpt_utils import EPS, IDATA, INFINITY, convertToCOO, convertToDense, extractRows, mapToCSC, scaleRows

# isort: off


class Optimizer(BaseSolver):
    def __init__(
        self,
        name: str,
        category: str,
        defaultOptions: Dict[str, Any] = {},
        informs: Dict[int, str] = {},
        options: Dict[str, Any] = {},
        checkDefaultOptions: bool = True,
        caseSensitiveOptions: bool = True,
        version: Optional[str] = None,
    ):
        """
        This is the base optimizer class that all optimizers inherit from.
        We define common methods here to avoid code duplication.

        Parameters
        ----------
        name : str
            Optimizer name
        category : str
            Typically local or global
        defaultOptions : dictionary
            A dictionary containing the default options
        informs : dict
            Dictionary of the inform codes
        """
        super().__init__(
            name,
            category,
            defaultOptions=defaultOptions,
            options=options,
            informs=informs,
            checkDefaultOptions=checkDefaultOptions,
            caseSensitiveOptions=caseSensitiveOptions,
        )
        # callCounter will be incremented after the function calls, iterCounters will be incremented before the calls.
        self.callCounter = 0  # counts all function calls (fobj, fcon, gobj, gcon)
        self.iterCounter = -1  # counts iteration(new x point)
        self.sens: Union[None, Callable, Gradient] = None
        self.optProb: Optimization
        self.version: Optional[str] = version

        # Default options:
        self.appendLinearConstraints: bool = False
        self.jacType: str = "dense"
        self.unconstrained: bool = False
        self.userObjTime: float = 0.0
        self.userSensTime: float = 0.0
        self.interfaceTime: float = 0.0
        self.userObjCalls: int = 0
        self.userSensCalls: int = 0
        self.storeSens: bool = True

        # Cache storage
        self.cache: Dict[str, Any] = {"x": None, "fobj": None, "fcon": None, "gobj": None, "gcon": None, "fail": None}

        # A second-level cache for optimizers that require callbacks
        # for each constraint. (eg. PSQP etc)
        self.storedData: Dict[str, Any] = {"x": None}

        # Store the Jacobian conversion maps
        self._jac_map_csr_to_csc = None

        # Initialize metadata
        self.metadata: Dict[str, Any] = {}
        self.startTime = None

    def _clearTimings(self):
        """Clear timings and call counters"""
        self.userObjTime = 0.0
        self.userSensTime = 0.0
        self.interfaceTime = 0.0
        self.userObjCalls = 0
        self.userSensCalls = 0

    def _setSens(self, sens: Union[None, str, Callable], sensStep: float, sensMode: str):
        """
        Common function to setup sens function
        """

        # If the sens parameter is None and the sens parameter in the
        # optProb is not None, use the optProb setting
        if sens is None and self.optProb.sens is not None:
            sens = self.optProb.sens

        # If we have SNOPT set derivative level to 3...it will be
        # reset if necessary
        if self.name in ["SNOPT"]:
            # SNOPT is the only one where None is ok.
            self.setOption("Derivative level", 3)

        # Next we determine what to what to do about
        # derivatives. We must have a function or we use FD or CS:
        if sens is None:
            if self.name in ["SNOPT"]:
                # SNOPT is the only one where None is ok.
                self.setOption("Derivative level", 0)
                self.sens = None
            else:
                raise ValueError(
                    "'None' value given for sens. "
                    + "Must be one of 'FD', 'FDR', 'CD', 'CDR', 'CS' or a user supplied function."
                )
        elif callable(sens):
            # We have function handle for gradients! Excellent!
            self.sens = sens
        elif sens.lower() in ["fd", "fdr", "cd", "cdr", "cs"]:
            # Create the gradient class that will operate just like if
            # the user supplied function
            self.sens = Gradient(self.optProb, sens.lower(), sensStep, sensMode, self.optProb.comm)
        else:
            raise ValueError(
                "Unknown value given for sens. Must be one of [None,'FD','FDR','CD','CDR','CS'] or a python function handle"
            )

    def _setHistory(self, storeHistory: str, hotStart: str):
        """
        Generic routine for setting up the hot start information

        Parameters
        ----------
        storeHistory : str
            File for possible history file. Or None if not writing file.

        hotStart : str
            Filename for history file for hot start
        """
        # we have to wrap the whole function
        # so it's parallel safe
        if self.optProb.comm.rank == 0:
            # By default no hot start
            self.hotStart = None

            # Determine if we want to do a hot start:
            if hotStart is not None:
                # Now, if if the hot start file and the history are
                # the SAME, we don't allow that. We will create a copy
                # of the hotStart file and use *that* instead.
                if storeHistory == hotStart:
                    if os.path.exists(hotStart):
                        fname = tempfile.mktemp()
                        shutil.copyfile(storeHistory, fname)
                        self.hotStart = History(fname, temp=True, flag="r")
                else:
                    if os.path.exists(hotStart):
                        self.hotStart = History(hotStart, temp=False, flag="r")
                    else:
                        pyOptSparseWarning("Hot start file does not exist. Performing a regular start")

            self.storeHistory = False
            if storeHistory:
                self.hist = History(storeHistory, flag="n", optProb=self.optProb)
                self.storeHistory = True

                if self.hotStart is not None:
                    # we set the DVs to the _initial_ values of the hotstart history file
                    # we need major=False here since not all optimizers support major iteration counting
                    # even though in theory the first call counter should always be major
                    init_DV = self.hotStart.getValues(
                        names=self.hotStart.getDVNames(), callCounters=[0], major=False, allowSens=True
                    )
                    self.optProb.setDVs(init_DV)
                    # we also save these metadata values
                    # into the new history file
                    for key in ["varInfo", "conInfo", "objInfo", "optProb"]:
                        val = self.hotStart.read(key)
                        if val is not None:
                            self.hist.writeData(key, val)
                    self._setMetadata()
                    self.hist.writeData("metadata", self.metadata)
        self.optProb.comm.Barrier()

    def _masterFunc(self, x: ndarray, evaluate: List[str]):
        """
        This is the master function that **ALL** optimizers call from
        the specific signature functions. The reason for this is that
        we can generically do the hot-start replay, history storage,
        timing and possibly caching once for all optimizers. It also
        takes care of the MPI communication that allows the optimizer
        to run on one process only, but within a larger MPI context.

        It does add one additional level of call, but we think it is
        well worth it for reduce code duplication

        Parameters
        ----------
        x : array
            This is the raw x-array data from the optimizer
        evaluate : list of strings
            This list contains at least one of 'fobj', 'fcon', 'gobj'
            or 'gcon'. This list tells this function which of the
            values is required on return
        """

        # Increment iteration counter if x is a new point
        if not np.isclose(x, self.cache["x"], atol=EPS, rtol=EPS).all():
            self.iterCounter += 1

        # We are hot starting, we should be able to read the required
        # information out of the hot start file, process it and then
        # fire it back to the specific optimizer
        timeA = time.time()
        if self.hotStart:
            # This is a very inexpensive check to see if point exists
            if self.hotStart.pointExists(self.callCounter):
                # Read the actual data for this point:
                data = self.hotStart.read(self.callCounter)

                # Get the x-value and (de)process
                xuser_ref = self.optProb.processXtoVec(data["xuser"])

                # Validated x-point point to use:
                xuser_vec = self.optProb._mapXtoUser(x)
                if np.isclose(xuser_vec, xuser_ref, rtol=EPS, atol=EPS).all():
                    # However, we may need a sens that *isn't* in the
                    # the dictionary:
                    funcs = None
                    funcsSens = None
                    validPoint = True
                    if "fobj" in evaluate or "fcon" in evaluate:
                        funcs = data["funcs"]

                    if "gobj" in evaluate or "gcon" in evaluate:
                        if "funcsSens" in data:
                            funcsSens = data["funcsSens"]
                        else:
                            validPoint = False

                    # Only continue if valid:
                    if validPoint:
                        if self.storeHistory:
                            # Just dump the (exact) dictionary back out:
                            data["isMajor"] = False
                            self.hist.write(self.callCounter, data)

                        fail = data["fail"]
                        returns = []

                        # Process constraints/objectives
                        if funcs is not None:
                            self.optProb.evaluateLinearConstraints(xuser_vec, funcs)
                            fcon = self.optProb.processContoVec(funcs)
                            fobj = self.optProb.processObjtoVec(funcs)
                            if "fobj" in evaluate:
                                returns.append(fobj)
                            if "fcon" in evaluate:
                                returns.append(fcon)

                        # Process gradients if we have them
                        if funcsSens is not None:
                            gobj = self.optProb.processObjectiveGradient(funcsSens)
                            gcon = self.optProb.processConstraintJacobian(funcsSens)
                            gcon = self._convertJacobian(gcon)

                            if "gobj" in evaluate:
                                returns.append(gobj)
                            if "gcon" in evaluate:
                                returns.append(gcon)

                        # Cache x because the iteration counter need this
                        self.cache["x"] = x.copy()

                        # We can now safely increment the call counter
                        self.callCounter += 1
                        returns.append(fail)
                        self.interfaceTime += time.time() - timeA
                        return returns
                    # end if (valid point -> all data present)
                # end if (x's match)
            # end if (point exists)

            # We have used up all the information in hot start so we
            # can close the hot start file
            self.hotStart.close()
            self.hotStart = None
        # end if (hot starting)

        # Now we have to actually run our function...this is where the
        # MPI gets a little tricky. Up until now, only the root proc
        # has called up to here...the rest of them are waiting at a
        # broadcast to know what to do.

        args = [x, evaluate]

        # Broadcast the type of call (0 means regular call)
        self.optProb.comm.bcast(0, root=0)

        # Now broadcast out the required arguments:
        self.optProb.comm.bcast(args)

        result = self._masterFunc2(*args)
        self.interfaceTime += time.time() - timeA
        return result

    def _masterFunc2(self, x, evaluate, writeHist=True):
        """
        Another shell function. This function is now actually called
        on all the processors.
        """

        # Our goal in this function is to return the values requested
        # in 'evaluate' for the corresponding x. We have to be a
        # little cheeky here since some optimizers will make multiple
        # call backs with the same x, one for the objective and one
        # for the constraint. We therefore at the end of each function
        # or sensitivity call we cache the x value and the fobj, fcon,
        # gobj, and gcon values such that on the next pass we can just
        # read them and return.

        xuser_vec = self.optProb._mapXtoUser(x)
        xuser = self.optProb.processXtoDict(xuser_vec)

        masterFail = 0

        # Set basic parameters in history
        hist = {"xuser": xuser}
        returns = []
        # Start with fobj:
        if "fobj" in evaluate:
            if not np.isclose(x, self.cache["x"], atol=EPS, rtol=EPS).all() or "funcs" not in self.cache:
                # The previous evaluated point is different than the point requested
                # OR this is a recursive call to _masterFunc2 from a gradient evaluation that occured
                # at the beginning of a hot started optimization
                timeA = time.time()
                args = self.optProb.objFun(xuser)
                if isinstance(args, tuple):
                    funcs = args[0]
                    fail = args[1]
                elif args is None:
                    raise ValueError(
                        "No return values from user supplied objective function. "
                        + "The function must return 'funcs' or 'funcs, fail'"
                    )
                else:
                    funcs = args
                    fail = 0

                self.userObjTime += time.time() - timeA
                self.userObjCalls += 1

                # Make sure the user-defined function does *not* return linear constraint values
                if self.callCounter == 0:
                    self._checkLinearConstraints(funcs)

                # Discard zero imaginary components in funcs
                for key, val in funcs.items():
                    funcs[key] = np.real(val)

                # Store user values
                self.cache["funcs"] = copy.deepcopy(funcs)

                # Process constraints/objectives
                self.optProb.evaluateLinearConstraints(xuser_vec, funcs)
                fcon = self.optProb.processContoVec(funcs)
                fobj = self.optProb.processObjtoVec(funcs)
                # Now clear out gobj and gcon in the cache since these
                # are out of date and set the current ones
                self.cache["gobj"] = None
                self.cache["gcon"] = None
                self.cache["x"] = x.copy()
                self.cache["fobj"] = copy.deepcopy(fobj)
                self.cache["fcon"] = copy.deepcopy(fcon)

                # Update fail flag
                masterFail = max(masterFail, fail)
                self.cache["fail"] = masterFail

            # fobj is now in cache
            returns.append(self.cache["fobj"])
            hist["funcs"] = self.cache["funcs"]

        if "fcon" in evaluate:
            if not np.isclose(x, self.cache["x"], atol=EPS, rtol=EPS).all() or "funcs" not in self.cache:
                # The previous evaluated point is different than the point requested
                # OR this is a recursive call to _masterFunc2 from a gradient evaluation that occured
                # at the beginning of a hot started optimization
                timeA = time.time()

                args = self.optProb.objFun(xuser)
                if isinstance(args, tuple):
                    funcs = args[0]
                    fail = args[1]
                elif args is None:
                    raise ValueError(
                        "No return values from user supplied objective function. "
                        + "The function must return 'funcs' *OR* 'funcs, fail'"
                    )
                else:
                    funcs = args
                    fail = 0

                self.userObjTime += time.time() - timeA
                self.userObjCalls += 1

                # Make sure the user-defined function does *not* return linear constraint values
                if self.callCounter == 0:
                    self._checkLinearConstraints(funcs)

                # Discard zero imaginary components in funcs
                for key, val in funcs.items():
                    funcs[key] = np.real(val)

                # Store user values
                self.cache["funcs"] = copy.deepcopy(funcs)

                # Process constraints/objectives
                self.optProb.evaluateLinearConstraints(xuser_vec, funcs)
                fcon = self.optProb.processContoVec(funcs)
                fobj = self.optProb.processObjtoVec(funcs)
                # Now clear out gobj and gcon in the cache since these
                # are out of date and set the current ones
                self.cache["gobj"] = None
                self.cache["gcon"] = None
                self.cache["x"] = x.copy()
                self.cache["fobj"] = copy.deepcopy(fobj)
                self.cache["fcon"] = copy.deepcopy(fcon)

                # Update fail flag
                masterFail = max(masterFail, fail)
                self.cache["fail"] = masterFail

            # fcon is now in cache
            returns.append(self.cache["fcon"])
            hist["funcs"] = self.cache["funcs"]

        if "gobj" in evaluate:
            if not np.isclose(x, self.cache["x"], atol=EPS, rtol=EPS).all() or "funcs" not in self.cache:
                # The previous evaluated point is different than the point requested for the derivative
                # OR this is the first call to _masterFunc2 in a hot started optimization
                # Recursively call the routine with ['fobj', 'fcon']
                _, _, fail = self._masterFunc2(x, ["fobj", "fcon"], writeHist=False)
                # We *don't* count that extra call, since that will
                # screw up the numbering...so we subtract the last call.
                self.callCounter -= 1
                # Update fail flag
                masterFail = max(masterFail, fail)
                self.cache["fail"] = masterFail
            # Now, the point has been evaluated correctly so we
            # determine if we have to run the sens calc:

            if self.cache["gobj"] is None:
                timeA = time.time()
                args = self.sens(xuser, self.cache["funcs"])

                if isinstance(args, tuple):
                    funcsSens = args[0]
                    fail = args[1]
                elif args is None:
                    raise ValueError(
                        "No return values from user supplied sensitivity function. "
                        + "The function must return 'funcsSens' or 'funcsSens, fail'"
                    )
                else:
                    funcsSens = args
                    fail = 0

                self.userSensTime += time.time() - timeA
                self.userSensCalls += 1

                # User values are stored immediately
                # deepcopy of the sens dictionary is slow, so just reference it
                # It shouldn't be modified until the next sensitivity call.
                self.cache["funcsSens"] = funcsSens

                # Process objective gradient for optimizer
                gobj = self.optProb.processObjectiveGradient(funcsSens)

                # Process constraint gradients for optimizer
                gcon = self.optProb.processConstraintJacobian(funcsSens)
                gcon = self._convertJacobian(gcon)

                # Set the cache values:
                self.cache["gobj"] = gobj.copy()
                self.cache["gcon"] = gcon.copy()

                # Update fail flag
                masterFail = max(masterFail, fail)
                self.cache["fail"] = masterFail

            # gobj is now in the cache
            returns.append(self.cache["gobj"])
            if self.storeSens:
                hist["funcsSens"] = self.cache["funcsSens"]

        if "gcon" in evaluate:
            if not np.isclose(x, self.cache["x"], atol=EPS, rtol=EPS).all() or "funcs" not in self.cache:
                # The previous evaluated point is different than the point requested for the derivative
                # OR this is the first call to _masterFunc2 in a hot started optimization
                # Recursively call the routine with ['fobj', 'fcon']
                _, _, fail = self._masterFunc2(x, ["fobj", "fcon"], writeHist=False)
                # We *don't* count that extra call, since that will
                # screw up the numbering...so we subtract the last call.
                self.callCounter -= 1
                # Update fail flag
                masterFail = max(masterFail, fail)
                self.cache["fail"] = masterFail
            # Now, the point has been evaluated correctly so we
            # determine if we have to run the sens calc:
            if self.cache["gcon"] is None:
                timeA = time.time()

                args = self.sens(xuser, self.cache["funcs"])

                if isinstance(args, tuple):
                    funcsSens = args[0]
                    fail = args[1]
                elif args is None:
                    raise ValueError(
                        "No return values from user supplied sensitivity function. "
                        + "The function must 'return 'funcsSens' or 'funcsSens, fail'"
                    )
                else:
                    funcsSens = args
                    fail = 0

                self.userSensTime += time.time() - timeA
                self.userSensCalls += 1

                # User values are stored immediately
                self.cache["funcsSens"] = funcsSens

                # Process objective gradient for optimizer
                gobj = self.optProb.processObjectiveGradient(funcsSens)

                # Process constraint gradients for optimizer
                gcon = self.optProb.processConstraintJacobian(funcsSens)
                gcon = self._convertJacobian(gcon)

                # Set cache values
                self.cache["gobj"] = gobj.copy()
                self.cache["gcon"] = gcon.copy()

                # Update fail flag
                masterFail = max(masterFail, fail)
                self.cache["fail"] = masterFail

            # gcon is now in the cache
            returns.append(self.cache["gcon"])
            if self.storeSens:
                hist["funcsSens"] = self.cache["funcsSens"]

        # Update the fail flag with any cached failure and put the fail flag in the history
        masterFail = max(self.cache["fail"], masterFail)
        hist["fail"] = masterFail

        # Put the iteration counter in the history
        hist["iter"] = self.iterCounter

        # timing
        hist["time"] = time.time() - self.startTime

        # Save information about major iteration counting (only matters for SNOPT).
        if self.name == "SNOPT":
            hist["isMajor"] = False  # this will be updated in _snstop if it is major
        else:
            hist["isMajor"] = True  # for other optimizers we assume everything's major

        # Add constraint and variable bounds at beginning of optimization.
        # This info is used for visualization using OptView.
        if self.callCounter == 0 and self.optProb.comm.rank == 0:
            conInfo = OrderedDict()
            varInfo = OrderedDict()
            objInfo = OrderedDict()

            # Cycle through constraints adding the bounds
            for key in self.optProb.constraints.keys():
                if not self.optProb.constraints[key].linear:
                    lower = self.optProb.constraints[key].lower
                    upper = self.optProb.constraints[key].upper
                    scale = self.optProb.constraints[key].scale
                    conInfo[key] = {"lower": lower, "upper": upper, "scale": scale}

            # Cycle through variables and add the bounds
            for dvGroup in self.optProb.variables:
                varInfo[dvGroup] = {"lower": [], "upper": [], "scale": []}
                for var in self.optProb.variables[dvGroup]:
                    if var.type == "c":
                        varInfo[dvGroup]["lower"].append(var.lower / var.scale)
                        varInfo[dvGroup]["upper"].append(var.upper / var.scale)
                        varInfo[dvGroup]["scale"].append(var.scale)

            for objKey in self.optProb.objectives.keys():
                objInfo[objKey] = {"scale": self.optProb.objectives[objKey].scale}

            # There is a special write for additional metadata
            if self.storeHistory:
                self.hist.writeData("varInfo", varInfo)
                self.hist.writeData("conInfo", conInfo)
                self.hist.writeData("objInfo", objInfo)
                self._setMetadata()
                self.hist.writeData("metadata", self.metadata)
                # we have to get rid of some callables in optProb before serialization
                optProb = copy.copy(self.optProb)
                optProb.objFun = None
                optProb.sens = None
                self.hist.writeData("optProb", optProb)

        # Write history if necessary
        if self.optProb.comm.rank == 0 and writeHist and self.storeHistory:
            self.hist.write(self.callCounter, hist)

        # We can now safely increment the call counter
        self.callCounter += 1

        # Tack the fail flag on at the end
        returns.append(masterFail)

        return returns

    def _internalEval(self, x):
        """
        Special internal evaluation for optimizers that have a
        separate callback for each constraint"""

        fobj, fcon, gobj, gcon, fail = self._masterFunc(x, ["fobj", "fcon", "gobj", "gcon"])

        self.storedData["x"] = x.copy()
        self.storedData["fobj"] = fobj
        self.storedData["fcon"] = fcon.copy()
        self.storedData["gobj"] = gobj.copy()
        self.storedData["gcon"] = gcon.copy()

    def _checkEval(self, x):
        """Special check to be used with _internalEval()"""
        if self.storedData["x"] is None:
            return True
        elif (self.storedData["x"] == x).all():
            return False
        else:
            return True

    def _convertJacobian(self, gcon_csr_in):
        """
        Convert gcon which is a coo matrix into the format we need.

        The returned Jacobian gcon is the data only, not a dictionary.
        """

        # Now, gcon is a CSR sparse matrix.  Depending on what the
        # optimizer wants, we will convert. The conceivable options
        # are: dense (most), csc (snopt), csr (???), or coo (IPOPT)

        if self.optProb.nCon > 0:
            # Extract the rows we need:
            gcon_csr = extractRows(gcon_csr_in, self.optProb.jacIndices)

            # Apply factor scaling because of constraint sign changes
            scaleRows(gcon_csr, self.optProb.fact)

            # Now convert to final format:
            if self.jacType == "dense2d":
                gcon = convertToDense(gcon_csr)
            elif self.jacType == "csc":
                if self._jac_map_csr_to_csc is None:
                    self._jac_map_csr_to_csc = mapToCSC(gcon_csr)
                gcon = gcon_csr["csr"][IDATA][self._jac_map_csr_to_csc[IDATA]]
            elif self.jacType == "csr":
                gcon = gcon_csr["csr"][IDATA]
            elif self.jacType == "coo":
                gcon = convertToCOO(gcon_csr)
                gcon = gcon["coo"][IDATA]
        if self.optProb.dummyConstraint:
            gcon = gcon_csr_in["csr"][IDATA]
        return gcon

    def _waitLoop(self):
        """Non-root processors go into this waiting loop while the
        root proc does all the work in the optimization algorithm"""

        mode = None
        info = None
        while True:
            # * Note*: No checks for MPI here since this code is
            # * only run in parallel, which assumes mpi4py is working

            # Receive mode and quit if mode is -1:
            mode = self.optProb.comm.bcast(mode, root=0)
            if mode == -1:
                break

            # Otherwise receive info from shell function
            info = self.optProb.comm.bcast(info, root=0)

            # Call the generic internal function. We don't care
            # about return values on these procs
            self._masterFunc2(*info)

    def _setInitialCacheValues(self):
        """
        Once we know that the optProb has been set, we populate the
        cache with a magic number. If the starting points for your
        optimization is -9999999999 then you out of luck!
        """
        self.cache["x"] = -999999999 * np.ones(self.optProb.ndvs)

    def _assembleContinuousVariables(self):
        """
        Utility function for assembling the design variables. Most
        optimizers here use continuous variables so this chunk of code
        can be reused.
        """
        blx = []
        bux = []
        xs = []
        for dvGroup in self.optProb.variables:
            for var in self.optProb.variables[dvGroup]:
                if var.type == "c":
                    blx.append(var.lower)
                    bux.append(var.upper)
                    xs.append(var.value)

                else:
                    raise ValueError(f"{self.name} cannot handle integer or discrete design variables")

        blx = np.array(blx)
        bux = np.array(bux)
        xs = np.array(xs)

        return blx, bux, xs

    def _assembleConstraints(self):
        """
        Utility function for assembling the design variables. Most
        optimizers here use continuous variables so this chunk of code
        can be reused.
        """

        # Constraints Handling -- make sure nonlinear constraints
        # go first -- this is particular to slsqp
        blc = []
        buc = []

        for key in self.optProb.constraints.keys():
            if not self.optProb.constraints[key].linear:
                blc.extend(self.optProb.constraints[key].lower)
                buc.extend(self.optProb.constraints[key].upper)

        for key in self.optProb.constraints.keys():
            if self.optProb.constraints[key].linear:
                blc.extend(self.optProb.constraints[key].lower)
                buc.extend(self.optProb.constraints[key].upper)

        if self.unconstrained:
            blc.append(-INFINITY)
            buc.append(INFINITY)

        ncon = len(blc)
        blc = np.array(blc)
        buc = np.array(buc)

        return ncon, blc, buc

    def _assembleObjective(self):
        """
        Utility function for assembling the objective, fetching the information in the Objective object within the Optimization class.
        Most optimizers use a single objective. In that case, the function will return a 0-length array (not a scalar).
        """

        nobj = len(self.optProb.objectives.keys())
        ff = []
        if nobj == 0:
            raise ValueError("No objective function was supplied! One can be added using a call to optProb.addObj()")
        for objKey in self.optProb.objectives:
            ff.append(self.optProb.objectives[objKey].value)

        return np.real(np.squeeze(ff))

    def _createSolution(self, optTime, sol_inform, obj, xopt, multipliers=None) -> Solution:
        """
        Generic routine to create the solution after an optimizer
        finishes.
        """
        fStar = self.optProb._mapObjtoUser(obj)
        # optionally convert to dict for multi-objective problems
        if isinstance(fStar, (list, np.ndarray)) and len(fStar) > 1:
            fStar = self.optProb.processObjtoDict(fStar, scaled=False)
        xuser = self.optProb._mapXtoUser(xopt)
        xStar = self.optProb.processXtoDict(xuser)

        if multipliers is not None:
            multipliers = self.optProb.processContoDict(multipliers, scaled=True, multipliers=True)

            # objective scaling
            if len(self.optProb.objectives.keys()) == 1:  # we assume there is only one objective
                obj = list(self.optProb.objectives.keys())[0]
                for con in multipliers.keys():
                    multipliers[con] /= self.optProb.objectives[obj].scale
        # construct info dict
        info = {
            "optTime": optTime,
            "userObjTime": self.userObjTime,
            "userSensTime": self.userSensTime,
            "userObjCalls": self.userObjCalls,
            "userSensCalls": self.userSensCalls,
            "interfaceTime": self.interfaceTime - self.userSensTime - self.userObjTime,
            "optCodeTime": optTime - self.interfaceTime,
        }
        sol = Solution(self.optProb, xStar, fStar, multipliers, sol_inform, info)

        return sol

    def _communicateSolution(self, sol: Optional[Solution]) -> Solution:
        """
        Broadcast the solution from the root proc back to everyone. We
        have to be a little careful since we can't in general
        broadcast the function and comm so we have to set manually after the broadcast.
        """

        if sol is not None:
            sol.comm = None
        commSol = self.optProb.comm.bcast(sol)
        commSol.objFun = self.optProb.objFun
        commSol.comm = self.optProb.comm

        return commSol

    def _setMetadata(self):
        """
        This function is used to set the self.metadata object.
        Importantly, this sets the startTime, so should be called just before the start
        of the optimization. endTime should be directly appended to the dictionary
        after optimization finishes.
        """
        options = copy.copy(self.options)
        # we remove entries which can't be stored properly in the history file
        if "snSTOP function handle" in options.keys():
            options.pop("snSTOP function handle")

        from .__init__ import __version__  # importing the pyoptsparse version

        # we store the metadata now, and write it later in optimizer calls
        # since we need the runtime at the end of optimization
        self.metadata = {
            "version": __version__,
            "optimizer": self.name,
            "optVersion": self.version,
            "optName": self.optProb.name,
            "nprocs": MPI.COMM_WORLD.size,
            "optOptions": options,
            "startTime": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        }

    def _on_setOption(self, name, value):
        """
        Set Optimizer Option Value (Optimizer Specific Routine)
        """
        pass

    def _checkLinearConstraints(self, funcs):
        """
        Makes sure that the user-defined obj/con function does not compute the linear constraint values
        because the linear constraints are exclusively defined by jac and bounds in addConGroup.
        """
        for conName in self.optProb.constraints:
            if self.optProb.constraints[conName].linear and conName in funcs:
                raise ValueError(
                    "Value for linear constraint returned from user obj function. Linear constraints "
                    + "are evaluated internally and should not be returned from the user's function."
                )

    def setOption(self, name, value=None):
        """
        Generic routine for all option setting. The routine does
        error checking on the type of the value.

        Parameters
        ----------
        name : str
            Name of the option to set
        value : varies
            Variable value to set.
        """

        super().setOption(name, value)
        # Now call the optimizer specific routine
        self._on_setOption(name, value)

    def _on_getOption(self, name):
        """
        Routine to be implemented by optimizer
        """
        pass

    def getOption(self, name):
        """
        Return the optimizer option value for name

        Parameters
        ----------
        name : str
            name of option for which to retrieve value

        Returns
        -------
        value : varies
            value of option for 'name'
        """

        # Call the optimizer specific routine
        self._on_getOption(name)

        return super().getOption(name)

    def _on_getInform(self, info):
        """
        Routine to be implemented by optimizer
        """
        try:
            return self.informs[info]
        except KeyError:
            return f"Unknown Exit Status, Exit Code {info}"

    def getInform(self, infocode=None):
        """
        Get optimizer result inform code at exit

        Parameters
        ----------
        infocode : int
            Integer information code
        """

        if infocode is None:
            return self.informs
        else:
            return self._on_getInform(infocode)


# =============================================================================
# Generic OPT Constructor
# =============================================================================

# List of optimizers as an enum
Optimizers = Enum("Optimizers", "SNOPT IPOPT SLSQP NLPQLP CONMIN NSGA2 PSQP ALPSO ParOpt")


def OPT(optName, *args, **kwargs):
    """
    This is a simple utility function that enables creating an
    optimizer based on the 'optName' string. This can be useful for
    doing optimization studies with respect to optimizer since you
    don't need massive if-statements.

    Parameters
    ----------
    optName : str or enum
       Either a string identifying the optimizer to create, e.g. "SNOPT", or
       an enum accessed via ``pyoptsparse.Optimizers``, e.g. ``Optimizers.SNOPT``.

    ``*args``, ``**kwargs`` : varies
       Passed to optimizer creation.

    Returns
    -------
    opt : pyOpt_optimizer inherited optimizer
       The desired optimizer
    """
    if isinstance(optName, str):
        optName = optName.lower()
    if optName == "snopt" or optName == Optimizers.SNOPT:
        from .pySNOPT.pySNOPT import SNOPT as opt
    elif optName == "ipopt" or optName == Optimizers.IPOPT:
        from .pyIPOPT.pyIPOPT import IPOPT as opt
    elif optName == "slsqp" or optName == Optimizers.SLSQP:
        from .pySLSQP.pySLSQP import SLSQP as opt
    elif optName == "nlpqlp" or optName == Optimizers.NLPQLP:
        from .pyNLPQLP.pyNLPQLP import NLPQLP as opt
    elif optName == "psqp" or optName == Optimizers.PSQP:
        from .pyPSQP.pyPSQP import PSQP as opt
    elif optName == "conmin" or optName == Optimizers.CONMIN:
        from .pyCONMIN.pyCONMIN import CONMIN as opt
    elif optName == "nsga2" or optName == Optimizers.NSGA2:
        from .pyNSGA2.pyNSGA2 import NSGA2 as opt
    elif optName == "alpso" or optName == Optimizers.ALPSO:
        from .pyALPSO.pyALPSO import ALPSO as opt
    elif optName == "paropt" or optName == Optimizers.ParOpt:
        from .pyParOpt.ParOpt import ParOpt as opt
    else:
        raise ValueError(
            (
                "The optimizer specified in 'optName' was not recognized. "
                + "The current list of supported optimizers is {}"
            ).format(list(map(str, Optimizers)))
        )

    # Create the optimizer and return it
    return opt(*args, **kwargs)


def list_optimizers() -> list[Optimizers]:
    """List all optimizers which were installed successfully and available for use"""
    all_optimizers = []
    for opt in Optimizers:
        try:
            OPT(opt)
            all_optimizers.append(opt)
        except ImportError:
            pass
    return all_optimizers
