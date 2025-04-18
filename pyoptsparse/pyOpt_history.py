# Standard Python modules
from collections import OrderedDict
import copy
import os

# External modules
import numpy as np
from sqlitedict import SqliteDict

# Local modules
from .pyOpt_error import pyOptSparseWarning
from .pyOpt_utils import EPS


class History:
    def __init__(self, fileName, optProb=None, temp=False, flag="r"):
        """
        This class is essentially a thin wrapper around a SqliteDict dictionary to facilitate
        operations with pyOptSparse

        Parameters
        ----------
        fileName : str
            File name for history file

        optProb : pyOpt_Optimization
            The optimization object

        temp : bool
            Flag to signify that the file should be deleted after it is
            closed

        flag : str
            String specifying the mode. Similar to what was used in shelve.
            ``n`` for a new database and ``r`` to read an existing one.
        """
        self.flag = flag
        if self.flag == "n":
            # If writing, we expliclty remove the file to
            # prevent old keys from "polluting" the new histrory
            if os.path.exists(fileName):
                os.remove(fileName)
            self.db = SqliteDict(fileName)
            self.optProb = optProb
        elif self.flag == "r":
            if os.path.exists(fileName):
                # we cast the db to OrderedDict so we do not have to
                # manually close the underlying db at the end
                self.db = OrderedDict(SqliteDict(fileName))
            else:
                raise FileNotFoundError(
                    f"The requested history file {fileName} to open in read-only mode does not exist."
                )
            self._processDB()
        else:
            raise ValueError("The flag argument to History must be 'r' or 'n'.")
        self.temp = temp
        self.fileName = fileName

    def close(self):
        """
        Close the underlying database.
        This should only be used in write mode. In read mode, we close the db
        during initialization.
        """
        if self.flag == "n":
            self.db.close()
            if self.temp:
                os.remove(self.fileName)

    def write(self, callCounter, data):
        """
        This is the main to write data. Basically, we just pass in
        the callCounter, the integer forming the key, and a dictionary
        which will be written to the key

        Parameters
        ----------
        callCounter : int
            the callCounter to write to
        data : dict
            the dictionary corresponding to the callCounter
        """

        # String key to database on disk
        key = str(callCounter)
        # if the point exists, we merely update with new data
        if self.pointExists(callCounter):
            oldData = self.read(callCounter)
            oldData.update(data)
            self.db[key] = oldData
        else:
            self.db[key] = data
        self.db["last"] = key
        self.db.sync()
        self.keys = list(self.db.keys())

    def writeData(self, key, data):
        """
        Write arbitrary `key:data` value to db.

        Parameters
        ----------
        key : str
            The key to be added to the history file
        data
            The data corresponding to the key. It can be anything as long as it is serializable
            in `sqlitedict`.
        """
        self.db[key] = data
        self.db.commit()
        self.keys = list(self.db.keys())

    def pointExists(self, callCounter):
        """
        Determine if callCounter is in the database

        Parameters
        ----------
        callCounter : int or str of int

        Returns
        -------
        bool
            True if the callCounter exists in the history file.
            False otherwise.
        """
        if isinstance(callCounter, int):
            callCounter = str(callCounter)
        return callCounter in self.keys

    def read(self, key):
        """
        Read data for an arbitrary key. Returns None if key is not found.
        If passing in a callCounter, it should be verified by calling pointExists() first.

        Parameters
        ----------
        key : str or int
            generic key[str] or callCounter[int]

        Returns
        -------
        dict
            The value corresponding to `key` is returned.
            If the key is not found, `None` is returned instead.
        """
        if isinstance(key, int):
            key = str(key)
        try:
            return self.db[key]
        except KeyError:
            return None

    def _searchCallCounter(self, x):
        """
        Searches through existing callCounters, and finds the one corresponding
        to an evaluation at the design vector `x`.
        returns `None` if the point did not match previous evaluations

        Parameters
        ----------
        x : ndarray
            The unscaled DV as a single array.

        Returns
        -------
        int
            The callCounter corresponding to the DV `x`.
            `None` is returned if no match was found.

        Notes
        -----
        The tolerance used for this is the value `numpy.finfo(numpy.float64).eps`.
        """
        last = int(self.db["last"])
        callCounter = None
        for i in range(last, 0, -1):
            key = str(i)
            xuser = self.optProb.processXtoVec(self.db[key]["xuser"])
            if np.isclose(xuser, x, atol=EPS, rtol=EPS).all() and "funcs" in self.db[key].keys():
                callCounter = i
                break
        return callCounter

    def _processDB(self):
        """
        Pre-processes the DB file and store various values into class attributes.
        These will be used later when calling self.getXX functions.
        """
        # Load any keys it happens to have:
        self.keys = list(self.db.keys())
        # load info
        self.DVInfo = self.read("varInfo")
        self.conInfo = self.read("conInfo")
        self.objInfo = self.read("objInfo")
        # metadata
        self.metadata = self.read("metadata")
        self.optProb = self.read("optProb")
        # load names
        self.DVNames = set(self.getDVNames())
        self.conNames = set(self.getConNames())
        self.objNames = set(self.getObjNames())

        # extract list of callCounters from self.keys
        # this just checks if each key contains only digits, then cast into int
        self.callCounters = sorted([x for x in self.keys if x.isdigit()], key=float)

        # extract all information stored in the call counters
        self.iterKeys = set()
        self.extraFuncsNames = set()
        for i in self.callCounters:
            val = self.read(i)
            self.iterKeys.update(val.keys())
            if "funcs" in val.keys():
                self.extraFuncsNames.update(val["funcs"].keys())
        # remove objective and constraint keys
        self.extraFuncsNames = self.extraFuncsNames.difference(self.conNames).difference(self.objNames)

        from .__init__ import __version__  # isort: skip

        if self.metadata["version"] != __version__:
            pyOptSparseWarning(
                "The version of pyoptsparse used to generate the history file does not match the one being run right now. There may be compatibility issues."
            )

    def getIterKeys(self):
        """
        Returns the keys available at each optimization iteration.
        This function is useful for inspecting the history file, to determine
        what information is saved at each iteration.

        Returns
        -------
        list of str
            A list containing the names of keys stored at each optimization iteration.
        """
        return copy.deepcopy(list(self.iterKeys))

    def getDVNames(self):
        """
        Returns the names of the DVs.

        Returns
        -------
        list of str
            A list containing the names of DVs.
        """
        # only do this if we open the file with 'r' flag
        if self.flag != "r":
            return
        return copy.deepcopy(list(self.DVInfo.keys()))

    def getConNames(self):
        """
        Returns the names of constraints.

        Returns
        -------
        list of str
            A list containing the names of constraints.
        """
        # only do this if we open the file with 'r' flag
        if self.flag != "r":
            return
        # we remove linear constraints
        conNames = [con for con in self.conInfo.keys() if not self.optProb.constraints[con].linear]
        return copy.deepcopy(conNames)

    def getObjNames(self):
        """
        Returns the names of the objectives.

        Returns
        -------
        list of str
            A list containing the names of objectives.

        Notes
        -----
        Recall that for the sake of generality, pyOptSparse allows for multiple objectives to be
        added. This feature is not used currently, but does make `ObjNames` follow the same structure
        as `ConNames` and `DVNames`.
        """
        # only do this if we open the file with 'r' flag
        if self.flag != "r":
            return
        return copy.deepcopy(list(self.objInfo.keys()))

    def getExtraFuncsNames(self):
        """
        Returns extra funcs names.
        These are extra key: value pairs stored in the ``funcs`` dictionary for each iteration, which are not used by the optimizer.

        Returns
        -------
        list of str
            A list containing the names of extra funcs keys.

        """
        # only do this if we open the file with 'r' flag
        if self.flag != "r":
            return
        return copy.deepcopy(list(self.extraFuncsNames))

    def getObjInfo(self, key=None):
        """
        Returns the `ObjInfo`, for all keys by default. A `key` parameter can also
        be supplied, to retrieve `ObjInfo` corresponding to specific keys.


        Parameters
        ----------
        key : str or list of str, optional
            Specifies for which obj to extract `ObjInfo`.

        Returns
        -------
        dict
            A dictionary containing ObjInfo. For a single key, the return is one level deeper.

        Notes
        -----
        Recall that for the sake of generality, pyOptSparse allows for multiple objectives to be
        added. This feature is not used currently, but does make `ObjInfo` follow the same structure
        as `ConInfo` and `DVInfo`.
        Because of this, it is recommended that this function be accessed using the optional `key`
        argument.
        """
        # only do this if we open the file with 'r' flag
        if self.flag != "r":
            return
        if key is not None:
            if isinstance(key, str):
                return copy.deepcopy(self.objInfo[key])
            elif isinstance(key, list):
                d = OrderedDict()
                for k in key:
                    d[k] = self.objInfo[k]
                return d
        else:
            return copy.deepcopy(self.objInfo)

    def getConInfo(self, key=None):
        """
        Returns the `ConInfo`, for all keys by default. A `key` parameter can also
        be supplied, to retrieve `ConInfo` corresponding to specific constraints.

        Parameters
        ----------
        key : str or list of str, optional
            Specifies for which constraint to extract `ConInfo`.

        Returns
        -------
        dict
            A dictionary containing ConInfo. For a single key, the return is one level deeper.
        """
        # only do this if we open the file with 'r' flag
        if self.flag != "r":
            return
        if key is not None:
            if isinstance(key, str):
                return copy.deepcopy(self.conInfo[key])
            elif isinstance(key, list):
                d = OrderedDict()
                for k in key:
                    d[k] = self.conInfo[k]
                return d
        else:
            return copy.deepcopy(self.conInfo)

    def getDVInfo(self, key=None):
        """
        Returns the `DVInfo`, for all keys by default. A `key` parameter can also
        be supplied, to retrieve `DVInfo` corresponding to specific DVs.

        Parameters
        ----------
        key : str or list of str, optional
            Specifies for which DV to extract `DVInfo`.

        Returns
        -------
        dict
            A dictionary containing DVInfo. For a single key, the return is one level deeper.
        """
        # only do this if we open the file with 'r' flag
        if self.flag != "r":
            return
        if key is not None:
            if isinstance(key, str):
                return copy.deepcopy(self.DVInfo[key])
            elif isinstance(key, list):
                d = OrderedDict()
                for k in key:
                    d[k] = self.DVInfo[k]
                return d
        else:
            return copy.deepcopy(self.DVInfo)

    def getMetadata(self):
        """
        Returns a copy of the metadata stored in the history file.

        Returns
        -------
        dict
            A dictionary containing the metadata.
        """
        # only do this if we open the file with 'r' flag
        if self.flag != "r":
            return
        return copy.deepcopy(self.metadata)

    def getOptProb(self):
        """
        Returns a copy of the optProb associated with the optimization.

        Returns
        -------
        optProb : pyOpt_optimization object
            The optProb associated with the optimization. This is taken from the history file,
            and therefore has the ``comm`` and ``objFun`` fields removed.
        """
        # only do this if we open the file with 'r' flag
        if self.flag != "r":
            return
        return copy.deepcopy(self.optProb)

    def _processIterDict(self, d, scale=False):
        """
        This function scales the value, where the factor is extracted from the
        `Info` dictionaries, according to "name"

        Parameters
        ----------
        d : dictionary
            The iteration dictionary, i.e. hist['0']
            This must be a function evaluation callCounter, and not a gradient callCounter.
        scale : bool
            Whether the returned values should be scaled.

        Returns
        -------
        conDict : dict
            A dictionary containing constraint values
        objDict : dict
            A dictionary containing objective values
        DVDict : dict
            A dictionary containing DV values

        These are all "flat" dictionaries, with simple key:value pairs.
        """
        conDict = {}
        objDict = {}
        # these require funcs which may not always be there
        if "funcs" in d.keys():
            for con in list(self.optProb.constraints.keys()):
                # linear constraints are not stored in funcs
                if not self.optProb.constraints[con].linear:
                    conDict[con] = d["funcs"][con]
                else:
                    # the linear constraints are removed from optProb so that scaling works
                    # without needing the linear constraints to be present
                    self.optProb.constraints.pop(con)

            for obj in self.objNames:
                objDict[obj] = d["funcs"][obj]

        # the DVs will always be there
        DVDict = {}
        for DV in self.DVNames:
            DVDict[DV] = d["xuser"][DV]
        if scale:
            conDict = self.optProb._mapContoOpt_Dict(conDict)
            objDict = self.optProb._mapObjtoOpt_Dict(objDict)
            DVDict = self.optProb._mapXtoOpt_Dict(DVDict)
        return conDict, objDict, DVDict

    def getCallCounters(self):
        """
        Returns a list of all call counters stored in the history file.

        Returns
        -------
        list
            a list of strings, each entry being a call counter.
        """
        # only do this if we open the file with 'r' flag
        if self.flag != "r":
            return
        return copy.deepcopy(self.callCounters)

    def getValues(self, names=None, callCounters=None, major=True, scale=False, stack=False, allowSens=False):
        """
        Parses an existing history file and returns a data dictionary used to post-process optimization results, containing the requested optimization iteration history.

        Parameters
        ----------
        names : list or str
            the values of interest, can be the name of any DV, objective or constraint,
            or a list of them. If None, all values are returned. This includes the DVs,
            funcs, and any values stored by the optimizer.

        callCounters : list of ints, can also contain 'last'
            a list of callCounters to extract information from.
            If the callCounter is invalid, i.e. outside the range or is a funcsSens evaluation, then it is skipped.
            'last' represents the last major iteration.
            If None, values from all callCounters are returned.

        major : bool
            flag to specify whether to include only major iterations.

        scale : bool
            flag to specify whether to apply scaling for the values. True means
            to return values that are scaled the same way as the actual optimization.

        stack : bool
            flag to specify whether the DV should be stacked into a single numpy array with
            the key `xuser`, or retain their separate DVGroups.

        allowSens: bool
            flag to specify whether gradient evaluation iterations are allowed.
            If true, it is up to the user to ensure that the callCounters specified contain the information requested.

        Returns
        -------
        dict
            a dictionary containing the information requested. The keys of the dictionary
            correspond to the `names` requested. Each value is a numpy array with the first
            dimension equal to the number of callCounters requested.

        Notes
        -----
        Regardless of the major flag, failed function evaluations are not returned.

        Examples
        --------
        First we can request DV history over all major iterations:

        >>> hist.getValues(names='xvars', major=True)
        {'xvars': array([[-2.00000000e+00,  1.00000000e+00],
            [-1.00000000e+00,  9.00000000e-01],
            [-5.00305827e-17,  4.21052632e-01],
            [ 1.73666171e-06,  4.21049838e-01],
            [ 9.08477459e-06,  4.21045261e-01],
            [ 5.00000000e-01,  2.84786405e-01],
            [ 5.00000000e-01,  5.57279939e-01],
            [ 5.00000000e-01,  2.00000000e+00]])}

        Next we can look at DV and optimality for the first and last iteration only:

        >>> hist.getValues(names=['xvars','optimality'],callCounters=[0,'last'])
        {'optimality': array([1.27895528, 0. ]),
        'xvars': array([[-2. , 1. ],
                        [ 0.5, 2. ]])}
        """
        # only do this if we open the file with 'r' flag
        if self.flag != "r":
            return

        allNames = (
            self.DVNames.union(self.conNames)
            .union(self.objNames)
            .union(self.iterKeys)
            .union(self.extraFuncsNames)
            .difference({"funcs", "funcsSens", "xuser"})
        )
        # cast string input into a single list
        if isinstance(names, str):
            names = {names}
        elif names is None:
            names = allNames
        else:
            names = set(names)
        if stack:
            allNames.add("xuser")
        # error if names isn't either a DV, con or obj
        if not names.issubset(allNames):
            raise KeyError(
                "The names provided are not one of DVNames, conNames or objNames.\n"
                + f"The names must be a subset of {allNames}"
            )
        DVsAsFuncs = self.DVNames.intersection(self.conNames)
        if len(DVsAsFuncs) > 0:
            ambiguousNames = names.intersection(DVsAsFuncs)
            if len(ambiguousNames) > 0:
                pyOptSparseWarning(
                    f"The names provided {ambiguousNames} is ambiguous, since it is both a DV as well as an objective/constraint. "
                    + "It is being assumed to be a DV. If it was set up via addDVsAsFunctions, then there's nothing to worry. "
                    + "Otherwise, consider renaming the variable or manually editing the history file."
                )

        if len(names.intersection(self.iterKeys)) > 0:
            if not major:
                pyOptSparseWarning(
                    "The major flag has been set to True, since some names specified only exist on major iterations."
                )
                major = True

        if stack:
            DVinNames = names.intersection(self.DVNames)
            for DV in DVinNames:
                names.remove(DV)
            names.add("xuser")
            pyOptSparseWarning(
                "The stack flag was set to True. Therefore all DV names have been removed, and replaced with a single key 'xuser'."
            )

        # set up dictionary to return
        data = {}
        # pre-allocate list for each input
        for name in names:
            data[name] = []

        # this flag is used for error printing only
        user_specified_callCounter = False
        if callCounters is not None:
            user_specified_callCounter = True
            if isinstance(callCounters, str):
                callCounters = [callCounters]
        else:
            callCounters = self.callCounters

        # parse the 'last' callCounter
        if "last" in callCounters:
            callCounters.append(self.read("last"))
            callCounters.remove("last")

        self._previousIterCounter = -1
        # loop over call counters, check if each counter is valid, and parse
        for i in callCounters:
            val = self._readValidCallCounter(i, user_specified_callCounter, allowSens, major)
            if val is not None:  # if i is valid
                conDict, objDict, DVDict = self._processIterDict(val, scale=scale)
                for name in names:
                    if name == "xuser":
                        data[name].append(self.optProb.processXtoVec(DVDict))
                    elif name in self.DVNames:
                        data[name].append(DVDict[name])
                    elif name in self.conNames:
                        data[name].append(conDict[name])
                    elif name in self.objNames:
                        data[name].append(objDict[name])
                    elif name in self.extraFuncsNames:
                        data[name].append(val["funcs"][name])
                    else:  # must be opt
                        data[name].append(val[name])

        # reshape lists into numpy arrays
        for name in names:
            # we just stack along axis 0
            if len(data[name]) > 0:
                data[name] = np.stack(data[name], axis=0)
            else:
                data[name] = np.array(data[name])
            # we cast 1D arrays to 2D, for scalar values
            if data[name].ndim == 1:
                data[name] = np.expand_dims(data[name], 1)

        # Raise warning for IPOPT's duplicated history
        if self.db["metadata"]["optimizer"] == "IPOPT" and "iter" not in self.db["0"].keys():
            pyOptSparseWarning(
                "The optimization history of IPOPT has duplicated entries at every iteration. "
                + "Fix the history manually, or re-run the optimization with a current version of pyOptSparse to generate a correct history file. "
            )
        return data

    def _readValidCallCounter(self, i, user_specified_callCounter, allowSens, major):
        """
        Checks whether a call counter is valid and read the data. The call counter is valid when it is
            1) inside the range of the history data,
            2) a function evaluation (i.e. not a sensitivity evaluation, except when `allowSens = True`),
            3) not a duplicated entry,
            4) not a failed function evaluation,
            5) a major iteration (only when `major = True`).

        Parameters
        ----------
        i : int
            call counter.

        user_specified_callCounter : bool
            flag to specify whether the call counter `i` is requested by a user or not.

        allowSens: bool
            flag to specify whether gradient evaluation iterations are allowed.

        major : bool
            flag to specify whether to include only major iterations.

        Returns
        -------
        val : dict or None
            information corresponding to the call counter `i`.
            If the call counter is not valid, `None` is returned instead.
        """

        if not self.pointExists(i):
            if user_specified_callCounter:
                # user specified a non-existent call counter
                pyOptSparseWarning(f"callCounter {i} was not found and is skipped!")
            return None
        else:
            val = self.read(i)

            # check if the callCounter is of a function call
            if not ("funcs" in val.keys() or allowSens):
                if user_specified_callCounter:
                    # user unintentionally specified a call counter for sensitivity
                    pyOptSparseWarning(
                        f"callCounter {i} did not contain a function evaluation and is skipped! "
                        + "Was it a gradient evaluation step?"
                    )
                return None
            else:
                # exclude the duplicated history (only when we have "iter" recorded)
                if "iter" in val.keys():
                    duplicate_flag = val["iter"] == self._previousIterCounter
                    self._previousIterCounter = val["iter"]  # update iterCounter for next i
                    if duplicate_flag and not user_specified_callCounter:
                        # this is a duplicate
                        return None
                # end if "iter" in val.keys()

                # check major/minor iteration, and if the call failed
                if ((major and val["isMajor"]) or not major) and not val["fail"]:
                    return val
                else:
                    return None
            # end if - ("funcs" in val.keys()
        # end if - pointExists

    def __del__(self):
        try:
            self.db.close()
            if self.temp:
                os.remove(self.fileName)
        except Exception:
            pass
