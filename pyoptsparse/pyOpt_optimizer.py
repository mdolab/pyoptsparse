#!/usr/bin/env python
# =============================================================================
# Imports
# =============================================================================
import os
import time
import copy
import numpy as np
from .pyOpt_gradient import Gradient
from .pyOpt_error import Error, pyOptSparseWarning
from .pyOpt_history import History
from .pyOpt_solution import Solution
from .pyOpt_optimization import INFINITY
from .pyOpt_utils import convertToDense, convertToCOO, extractRows, mapToCSC, scaleRows, IDATA
from collections import OrderedDict
import datetime
from .pyOpt_MPI import MPI

eps = np.finfo(np.float64).eps

# =============================================================================
# Optimizer Class
# =============================================================================
class Optimizer(object):
    def __init__(self, name=None, category=None, defOptions=None, informs=None, **kwargs):
        """
        This is the base optimizer class that all optimizers inherit from.
        We define common methods here to avoid code duplication.

        Parameters
        ----------
        name : str
            Optimizer name
        category : str
            Typically local or global
        defOptions : dictionary
            A dictionary containing the default options
        informs : dict
            Dictionary of the inform codes
        """
        self.name = name
        self.category = category
        self.options = {}
        self.options["defaults"] = defOptions
        self.informs = informs
        self.callCounter = 0
        self.sens = None
        # Initialize Options
        for key in defOptions:
            self.options[key] = defOptions[key]

        koptions = kwargs.pop("options", {})
        for key in koptions:
            self.setOption(key, koptions[key])

        self.optProb = None
        # Default options:
        self.appendLinearConstraints = False
        self.jacType = "dense"
        self.unconstrained = False
        self.userObjTime = 0.0
        self.userSensTime = 0.0
        self.interfaceTime = 0.0
        self.userObjCalls = 0
        self.userSensCalls = 0
        self.storeSens = True

        # Cache storage
        self.cache = {"x": None, "fobj": None, "fcon": None, "gobj": None, "gcon": None}

        # A second-level cache for optimizers that require callbacks
        # for each constraint. (eg. PSQP etc)
        self.storedData = {}
        self.storedData["x"] = None

        # Store the Jacobian conversion maps
        self._jac_map_csr_to_csc = None

    def _clearTimings(self):
        """Clear timings and call counters"""
        self.userObjTime = 0.0
        self.userSensTime = 0.0
        self.interfaceTime = 0.0
        self.userObjCalls = 0
        self.userSensCalls = 0

    def _setSens(self, sens, sensStep, sensMode):
        """
        Common function to setup sens function
        """

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
                raise Error(
                    (
                        "'None' value given for sens. "
                        + "Must be one of 'FD', 'FDR', 'CD', 'CDR', 'CS' or a user supplied function."
                    )
                )
        elif hasattr(sens, "__call__"):
            # We have function handle for gradients! Excellent!
            self.sens = sens
        elif sens.lower() in ["fd", "fdr", "cd", "cdr", "cs"]:
            # Create the gradient class that will operate just like if
            # the user supplied function
            self.sens = Gradient(self.optProb, sens.lower(), sensStep, sensMode, self.optProb.comm)
        else:
            raise Error(
                "Unknown value given for sens. Must be one of [None,'FD','FDR','CD','CDR','CS'] or a python function handle"
            )

    def _setHistory(self, storeHistory, hotStart):
        """
        Generic routine for setting up the hot start information

        Parameters
        ----------
        storeHistory : str
            File for possible history file. Or None if not writing file.

        hotStart : str
            Filename for history file for hot start
        """
        # By default no hot start
        self.hotStart = None

        # Determine if we want to do a hot start:
        if hotStart is not None:
            # Now, if if the hot start file and the history are
            # the SAME, we don't allow that. We will create a copy
            # of the hotStart file and use *that* instead.
            import tempfile
            import shutil

            if storeHistory == hotStart:
                if os.path.exists(hotStart):
                    fname = tempfile.mktemp()
                    shutil.copyfile(storeHistory, fname)
                    self.hotStart = History(fname, temp=True, flag="r")
            else:
                if os.path.exists(hotStart):
                    self.hotStart = History(hotStart, temp=False, flag="r")
                else:
                    pyOptSparseWarning(
                        "Hot start file does not exist. \
                    Performing a regular start"
                    )

        self.storeHistory = False
        if storeHistory:
            self.hist = History(storeHistory, flag="n", optProb=self.optProb)
            self.storeHistory = True

            if hotStart is not None:
                for key in ["varInfo", "conInfo", "objInfo", "optProb"]:
                    val = self.hotStart.read(key)
                    if val is not None:
                        self.hist.writeData(key, val)
                self._setMetadata()
                self.hist.writeData("metadata", self.metadata)

    def _masterFunc(self, x, evaluate):
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
                if np.isclose(xuser_vec, xuser_ref, rtol=eps, atol=eps).all():

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
            if not np.isclose(x, self.cache["x"], atol=eps, rtol=eps).all():
                timeA = time.time()
                args = self.optProb.objFun(xuser)
                if isinstance(args, tuple):
                    funcs = args[0]
                    fail = args[1]
                elif args is None:
                    raise Error(
                        (
                            "No return values from user supplied objective function. "
                            + "The function must return 'funcs' or 'funcs, fail'"
                        )
                    )
                else:
                    funcs = args
                    fail = 0

                self.userObjTime += time.time() - timeA
                self.userObjCalls += 1
                # User values stored is immediately
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

            # fobj is now in cache
            returns.append(self.cache["fobj"])
            hist["funcs"] = self.cache["funcs"]

        if "fcon" in evaluate:
            if not np.isclose(x, self.cache["x"], atol=eps, rtol=eps).all():
                timeA = time.time()

                args = self.optProb.objFun(xuser)
                if isinstance(args, tuple):
                    funcs = args[0]
                    fail = args[1]
                elif args is None:
                    raise Error(
                        (
                            "No return values from user supplied objective function. "
                            + "The function must return 'funcs' *OR* 'funcs, fail'"
                        )
                    )
                else:
                    funcs = args
                    fail = 0

                self.userObjTime += time.time() - timeA
                self.userObjCalls += 1
                # User values stored is immediately
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

            # fcon is now in cache
            returns.append(self.cache["fcon"])
            hist["funcs"] = self.cache["funcs"]

        if "gobj" in evaluate:
            if not np.isclose(x, self.cache["x"], atol=eps, rtol=eps).all():
                # Previous evaluated point is *different* than the
                # point requested for the derivative. Recursively call
                # the routine with ['fobj', and 'fcon']
                self._masterFunc2(x, ["fobj", "fcon"], writeHist=False)
                # We *don't* count that extra call, since that will
                # screw up the numbering...so we subtract the last call.
                self.callCounter -= 1
            # Now, the point has been evaluated correctly so we
            # determine if we have to run the sens calc:

            if self.cache["gobj"] is None:
                timeA = time.time()
                args = self.sens(xuser, self.cache["funcs"])

                if isinstance(args, tuple):
                    funcsSens = args[0]
                    fail = args[1]
                elif args is None:
                    raise Error(
                        (
                            "No return values from user supplied sensitivity function. "
                            + "The function must return 'funcsSens' or 'funcsSens, fail'"
                        )
                    )
                else:
                    funcsSens = args
                    fail = 0

                self.userSensTime += time.time() - timeA
                self.userSensCalls += 1

                # User values are stored is immediately
                self.cache["funcsSens"] = copy.deepcopy(funcsSens)

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

            # gobj is now in the cache
            returns.append(self.cache["gobj"])
            if self.storeSens:
                hist["funcsSens"] = self.cache["funcsSens"]

        if "gcon" in evaluate:
            if not np.isclose(x, self.cache["x"], atol=eps, rtol=eps).all():
                # Previous evaluated point is *different* than the
                # point requested for the derivative. Recursively call
                # the routine with ['fobj', and 'fcon']
                self._masterFunc2(x, ["fobj", "fcon"], writeHist=False)
                # We *don't* count that extra call, since that will
                # screw up the numbering...so we subtract the last call.
                self.callCounter -= 1
            # Now, the point has been evaluated correctly so we
            # determine if we have to run the sens calc:
            if self.cache["gcon"] is None:
                timeA = time.time()

                args = self.sens(xuser, self.cache["funcs"])

                if isinstance(args, tuple):
                    funcsSens = args[0]
                    fail = args[1]
                elif args is None:
                    raise Error(
                        (
                            "No return values from user supplied sensitivity function. "
                            + "The function must 'return 'funcsSens' or 'funcsSens, fail'"
                        )
                    )
                else:
                    funcsSens = args
                    fail = 0

                self.userSensTime += time.time() - timeA
                self.userSensCalls += 1
                # User values stored is immediately
                self.cache["funcsSens"] = copy.deepcopy(funcsSens)

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

            # gcon is now in the cache
            returns.append(self.cache["gcon"])
            if self.storeSens:
                hist["funcsSens"] = self.cache["funcsSens"]

        # Put the fail flag in the history:
        hist["fail"] = masterFail

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

            # There is a special write for the bounds data
            if self.storeHistory:
                self.hist.writeData("varInfo", varInfo)
                self.hist.writeData("conInfo", conInfo)
                self.hist.writeData("objInfo", objInfo)
                self._setMetadata()
                self.hist.writeData("metadata", self.metadata)
                self.hist.writeData("optProb", self.optProb)

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
                pass
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
                    raise Error("%s cannot handle integer or discrete design variables" % self.name)

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
        Utility function for assembling the design variables. Most
        optimizers here use continuous variables so this chunk of code
        can be reused.
        """

        nobj = len(self.optProb.objectives.keys())
        ff = []
        if nobj == 0:
            raise Error("No objective function was supplied! One can be added using a call to optProb.addObj()")
        for objKey in self.optProb.objectives:
            ff.append(self.optProb.objectives[objKey].value)

        return np.real(np.squeeze(ff))

    def _createSolution(self, optTime, sol_inform, obj, xopt, multipliers=None):
        """
        Generic routine to create the solution after an optimizer
        finishes.
        """
        sol = Solution(self.optProb, optTime, sol_inform)
        sol.userObjTime = self.userObjTime
        sol.userSensTime = self.userSensTime
        sol.userObjCalls = self.userObjCalls
        sol.userSensCalls = self.userSensCalls
        sol.interfaceTime = self.interfaceTime - self.userSensTime - self.userObjTime
        sol.optCodeTime = sol.optTime - self.interfaceTime
        sol.fStar = obj  # FIXME: this doesn't work, at least for SNOPT
        xuser = self.optProb._mapXtoUser(xopt)
        sol.xStar = self.optProb.processXtoDict(xuser)

        # Now set the x-values:
        i = 0
        for dvGroup in sol.variables:
            for var in sol.variables[dvGroup]:
                var.value = xopt[i]
                i += 1

        if multipliers is not None:
            multipliers = self.optProb.processContoDict(multipliers, scaled=True, multipliers=True)

            # objective scaling
            if len(self.optProb.objectives.keys()) == 1:  # we assume there is only one objective
                obj = list(self.optProb.objectives.keys())[0]
                for con in multipliers.keys():
                    multipliers[con] /= self.optProb.objectives[obj].scale
            sol.lambdaStar = multipliers
        return sol

    def _communicateSolution(self, sol):
        """
        Broadcast the solution from the root proc back to everyone. We
        have to be a little careful since we can't in general
        broadcast the function and comm so we have to set manually after the broadcast.
        """

        if sol is not None:
            sol.comm = None
        sol = self.optProb.comm.bcast(sol)
        sol.objFun = self.optProb.objFun
        sol.comm = self.optProb.comm

        return sol

    def _setMetadata(self):
        """
        This function is used to set the self.metadata object.
        Importantly, this sets the startTime, so should be called just before the start
        of the optimization. endTime should be directly appended to the dictionary
        after optimization finishes.
        """
        options = copy.deepcopy(self.options)
        options.pop("defaults")  # remove the default list
        # we retrieve only the second item which is the actual value
        for key, val in options.items():
            options[key] = val[1]

        from .__init__ import __version__  # importing the pyoptsparse version

        # we store the metadata now, and write it later in optimizer calls
        # since we need the runtime at the end of optimization
        self.metadata = {
            "version": __version__,
            "optimizer": self.name,
            "optName": self.optProb.name,
            "nprocs": MPI.COMM_WORLD.size,
            "optOptions": options,
            "startTime": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        }

    def _on_setOption(self, name, value):
        """
        Set Optimizer Option Value (Optimizer Specific Routine)
        """
        raise Error("This optimizer has not implemented _on_setOption")

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

        if name in self.options["defaults"]:
            if type(value) == self.options["defaults"][name][0]:
                self.options[name] = [type(value), value]
            else:
                raise Error(
                    "Value type for option {} was incorrect. It was expecting type '{}' by received type '{}'".format(
                        name, self.options["defaults"][name][0], type(value)
                    )
                )
        else:
            raise Error("Received an unknown option: %s" % repr(name))

        # Now call the optimizer specific routine
        self._on_setOption(name, value)

    def _on_getOption(self, name):
        """
        Routine to be implemented by optimizer
        """
        raise Error("This optimizer has not implemented _on_getOption")

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

        if name in self.options["defaults"]:
            return self.options[name][1]
        else:
            raise Error("Received an unknown option: %s." % repr(name))

        # Now call the optimizer specific routine
        self._on_getOption(name)

    def _on_getInform(self, info):
        """
        Routine to be implemented by optimizer
        """
        raise Error("This optimizer has not implemented _on_getInform")

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


def OPT(optName, *args, **kwargs):
    """
    This is a simple utility function that enables creating an
    optimizer based on the 'optName' string. This can be useful for
    doing optimization studies with respect to optimizer since you
    don't need massive if-statements.

    Parameters
    ----------
    optName : str
       String identifying the optimizer to create

    *args, **kwargs : varies
       Passed to optimizer creation.

    Returns
    -------
    opt : pyOpt_optimizer inherited optimizer
       The desired optimizer
       """

    optName = optName.lower()
    optList = ["snopt", "ipopt", "slsqp", "nlpqlp", "conmin", "nsga2", "psqp", "alpso", "paropt"]
    if optName == "snopt":
        from .pySNOPT.pySNOPT import SNOPT as opt
    elif optName == "ipopt":
        from .pyIPOPT.pyIPOPT import IPOPT as opt
    elif optName == "slsqp":
        from .pySLSQP.pySLSQP import SLSQP as opt
    elif optName == "nlpqlp":
        from .pyNLPQLP.pyNLPQLP import NLPQLP as opt
    elif optName == "psqp":
        from .pyPSQP.pyPSQP import PSQP as opt
    elif optName == "conmin":
        from .pyCONMIN.pyCONMIN import CONMIN as opt
    elif optName == "nsga2":
        from .pyNSGA2.pyNSGA2 import NSGA2 as opt
    elif optName == "alpso":
        from .pyALPSO.pyALPSO import ALPSO as opt
    # elif optName == 'nomad':
    #     from .pyNOMAD.pyNOMAD import NOMAD as opt
    elif optName == "paropt":
        from .pyParOpt.ParOpt import ParOpt as opt
    else:
        raise Error(
            "The optimizer specified in 'optName' was \
not recognized. The current list of supported optimizers is: %s"
            % repr(optList)
        )

    # Create the optimizer and return it
    return opt(*args, **kwargs)
