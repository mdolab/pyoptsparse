#!/usr/bin/env python
# =============================================================================
# External Python modules
# =============================================================================
import os
import copy
import numpy as np
from .pyOpt_error import Error, pyOptSparseWarning
from sqlitedict import SqliteDict
from collections import OrderedDict

eps = np.finfo(np.float64).eps
# =============================================================================
# History Class
# =============================================================================
class History(object):
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
                raise Error("The requested history file %s to open in read-only mode does not exist." % fileName)
            self._processDB()
        else:
            raise Error("The flag argument to History must be 'r' or 'n'.")
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
        key = "%d" % callCounter
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
            key = "%d" % i
            xuser = self.optProb.processXtoVec(self.db[key]["xuser"])
            if np.isclose(xuser, x, atol=eps, rtol=eps).all() and "funcs" in self.db[key].keys():
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
        for i in self.callCounters:
            val = self.read(i)
            self.iterKeys.update(val.keys())

        from .__init__ import __version__

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
        conNames = list(self.conInfo.keys())
        for con in conNames:
            if self.optProb.constraints[con].linear:
                conNames.remove(con)
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
        for con in self.conNames:
            # linear constraints are not stored in funcs
            if not self.optProb.constraints[con].linear:
                conDict[con] = d["funcs"][con]
        objDict = {}
        for obj in self.objNames:
            objDict[obj] = d["funcs"][obj]
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

    def getValues(self, names=None, callCounters=None, major=True, scale=False, stack=False):
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
            .difference(set(["funcs", "funcsSens", "xuser"]))
        )
        # cast string input into a single list
        if isinstance(names, str):
            names = set([names])
        elif names is None:
            names = allNames
        else:
            names = set(names)
        if stack:
            allNames.add("xuser")
        # error if names isn't either a DV, con or obj
        if not names.issubset(allNames):
            raise Error(
                "The names provided are not one of DVNames, conNames or objNames.\n\
                The names must be a subset of {}".format(
                    allNames
                )
            )
        DVsAsFuncs = self.DVNames.intersection(self.conNames)
        if len(DVsAsFuncs) > 0:
            ambiguousNames = names.intersection(DVsAsFuncs)
            if len(ambiguousNames) > 0:
                pyOptSparseWarning(
                    "The names provided {} is ambiguous, since it is both a DV as well as an objective/constraint. It is being assumed to be a DV. If it was set up via addDVsAsFunctions, then there's nothing to worry. Otherwise, consider renaming the variable or manually editing the history file.".format(
                        ambiguousNames
                    )
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

        for i in callCounters:
            if self.pointExists(i):
                val = self.read(i)
                if "funcs" in val.keys():  # we have function evaluation
                    if ((major and val["isMajor"]) or not major) and not val["fail"]:
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
                            else:  # must be opt
                                data[name].append(val[name])
                    elif val["fail"] and user_specified_callCounter:
                        pyOptSparseWarning(
                            ("callCounter {} contained a failed function evaluation and is skipped!").format(i)
                        )
                elif user_specified_callCounter:
                    pyOptSparseWarning(
                        (
                            "callCounter {} did not contain a function evaluation and is skipped! Was it a gradient evaluation step?"
                        ).format(i)
                    )
            elif user_specified_callCounter:
                pyOptSparseWarning(("callCounter {} was not found and is skipped!").format(i))
        # reshape lists into numpy arrays
        for name in names:
            # we just stack along axis 0
            data[name] = np.stack(data[name], axis=0)
            # we cast 1D arrays to 2D, for scalar values
            if data[name].ndim == 1:
                data[name] = np.expand_dims(data[name], 1)

        return data

    def __del__(self):
        try:
            self.db.close()
            if self.temp:
                os.remove(self.fileName)
        except:  # noqa: E722
            pass


# ==============================================================================
# Optimizer History Test
# ==============================================================================
if __name__ == "__main__":

    # Test Optimizer History
    print("Testing Optimizer History...")
