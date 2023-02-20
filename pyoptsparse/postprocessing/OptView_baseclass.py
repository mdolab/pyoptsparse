"""
Shared base class for both OptView and OptView_dash.
This reduces code duplication by having both OptViews read from this baseclass.

John Jasa 2015-2019

"""

# Standard Python modules
import shelve

# External modules
import numpy as np
from sqlitedict import SqliteDict

# Local modules
from ..pyOpt_error import pyOptSparseWarning


class OVBaseClass:

    """
    Container for display parameters, properties, and objects.
    This includes a canvas for MPL plots and a bottom area with widgets.
    """

    def OptimizationHistory(self):
        """
        Reads in database history file and stores contents.
        Function information is stored as a dict in func_data,
        variable information is stored as a dict in var_data,
        and bounds information is stored as a dict in bounds.
        """

        # Initialize dictionaries for design variables and unknowns.
        # The data is saved redundantly in dicts for all iterations and then
        # for major iterations as well.
        self.func_data_all = {}
        self.func_data_major = {}
        self.var_data_all = {}
        self.var_data_major = {}
        db = {}
        self.num_iter = 0

        # Loop over each history file name provided by the user.
        for histIndex, histFileName in enumerate(self.histList):
            # If they only have one history file, we don't change the keys' names
            if len(self.histList) == 1:
                histIndex = ""
            else:  # If multiple history files, append letters to the keys,
                # such that 'key' becomes 'key_A', 'key_B', etc
                histIndex = "_" + chr(histIndex + ord("A"))
            self.histIndex = histIndex

            try:  # This is the classic method of storing history files
                db = shelve.open(histFileName, "r")
                OpenMDAO = False
            except Exception:  # Bare except because error is not in standard Python.
                # If the db has the 'iterations' tag, it's an OpenMDAO db.
                db = SqliteDict(histFileName, "iterations")
                OpenMDAO = True

                # Need to do this since in py3 db.keys() is a generator object
                keys = [i for i in db.keys()]

                # If it has no 'iterations' tag, it's a pyOptSparse db.
                if keys == []:
                    OpenMDAO = False
                    db = SqliteDict(histFileName)

            # Specific instructions for OpenMDAO databases
            if OpenMDAO:
                # Get the number of iterations by looking at the largest number
                # in the split string names for each entry in the db
                for string in db.keys():
                    string = string.split("|")

                nkey = int(string[-1])
                self.solver_name = string[0]

                # Initalize a list detailing if the iterations are major or minor
                self.iter_type = np.zeros(nkey)

                # Get the keys of the database where derivatives were evaluated.
                # These correspond to major iterations, while no derivative
                # info is calculated for gradient-free linesearches.
                deriv_keys = SqliteDict(histFileName, "derivs").keys()
                self.deriv_keys = [int(key.split("|")[-1]) for key in deriv_keys]

                # Save information from the history file for the funcs.
                self.DetermineMajorIterations(db, OpenMDAO=OpenMDAO)

                # Save information from the history file for the unknowns.
                self.SaveDBData(db, self.func_data_all, self.func_data_major, OpenMDAO=OpenMDAO, data_str="Unknowns")

                # Save information from the history file for the design variables.
                self.SaveDBData(db, self.var_data_all, self.var_data_major, OpenMDAO=OpenMDAO, data_str="Parameters")

                # Add labels to OpenMDAO variables.
                # Corresponds to constraints, design variables, and objective.
                try:
                    db = SqliteDict(histFileName, "metadata")
                    self.SaveOpenMDAOData(db)

                except KeyError:  # Skip metadata info if not included in OpenMDAO hist file
                    pass

            else:
                # Get the number of iterations
                nkey = int(db["last"]) + 1
                self.nkey = nkey

                # Initalize a list detailing if the iterations are major or minor
                # 1 = major, 2 = minor, 0 = sensitivity (or duplicated info by IPOPT)
                # The entries whose iter_type = 0 will be ignored.
                self.iter_type = np.zeros(nkey)

                # Check to see if there is bounds information in the db file.
                # If so, add them to self.bounds to plot later.
                try:
                    try:
                        info_dict = db["varInfo"].copy()
                        info_dict.update(db["conInfo"])
                        scale_info = True
                    except KeyError:
                        self.warning_display(
                            "This is an older optimization history file.\n"
                            + "Only bounds information has been stored, not scalar info."
                        )
                        info_dict = db["varBounds"].copy()
                        info_dict.update(db["conBounds"])
                        scale_info = False

                    # Got to be a little tricky here since we're modifying
                    # info_dict; if we simply loop over it with the generator
                    # from Python3, it will contain the new keys and then the
                    # names will be mangled incorrectly.
                    bounds_dict = {}
                    scaling_dict = {}
                    for key in info_dict.keys():
                        bounds_dict[key + histIndex] = {
                            "lower": info_dict[key]["lower"],
                            "upper": info_dict[key]["upper"],
                        }
                        if scale_info:
                            scaling_dict[key + histIndex] = info_dict[key]["scale"]

                    self.bounds.update(bounds_dict)
                    if scale_info:
                        self.scaling.update(scaling_dict)
                except KeyError:
                    pass

                # Check to see if there is proper saved info about iter type
                if "isMajor" in db["0"].keys():
                    self.storedIters = True
                else:
                    self.storedIters = False

                # Raise warning for IPOPT's duplicated history
                if "metadata" in db and db["metadata"]["optimizer"] == "IPOPT" and "iter" not in db["0"].keys():
                    pyOptSparseWarning(
                        "The optimization history file has duplicated entries at every iteration, and the OptView plot is not correct. "
                        + "Re-run the optimization with a current version of pyOptSparse to generate a correct history file."
                    )

                # Save information from the history file for the funcs.
                self.DetermineMajorIterations(db, OpenMDAO=OpenMDAO)

                # Save information from the history file for the funcs.
                self.SaveDBData(db, self.func_data_all, self.func_data_major, OpenMDAO=OpenMDAO, data_str="funcs")

                # Save information from the history file for the design variables.
                self.SaveDBData(db, self.var_data_all, self.var_data_major, OpenMDAO=OpenMDAO, data_str="xuser")

        # Set the initial dictionaries to reference all iterations.
        # Later this can be set to reference only the major iterations.
        self.func_data = self.func_data_all
        self.var_data = self.var_data_all

        # Find the maximum length of any variable in the dictionaries and
        # save this as the number of iterations.
        for data_dict in [self.func_data, self.var_data]:
            for key in data_dict.keys():
                length = len(data_dict[key])
                if length > self.num_iter:
                    self.num_iter = length

    def DetermineMajorIterations(self, db, OpenMDAO):
        if not OpenMDAO:
            previousIterCounter = -1

            # Loop over each optimization call
            for i, iter_type in enumerate(self.iter_type):
                # If this is an OpenMDAO file, the keys are of the format
                # 'rank0:SNOPT|1', etc
                key = "%d" % i

                # Only actual optimization iterations have 'funcs' in them.
                # pyOptSparse saves info for two iterations for every
                # actual major iteration. In particular, one has funcs
                # and the next has funcsSens, but they're both part of the
                # same major iteration.
                # For IPOPT, it saves info for four calls for every
                # actual major iteration: objective, constraints,
                # and sensitivities of each.

                if "funcs" in db[key].keys():
                    # check if this entry is duplicated info. Only relevant for IPOPT.
                    # Note: old hist files don't have "iter"
                    if "iter" in db[key].keys() and db[key]["iter"] == previousIterCounter:
                        # duplicated info
                        self.iter_type[i] = 0

                    # if we did not store major iteration info, everything's major
                    elif not self.storedIters:
                        self.iter_type[i] = 1
                    # this is major iteration
                    elif self.storedIters and db[key]["isMajor"]:
                        self.iter_type[i] = 1
                    else:
                        self.iter_type[i] = 2

                    if "iter" in db[key].keys():
                        previousIterCounter = db[key]["iter"]
                else:
                    self.iter_type[i] = 0  # this is not a real iteration,
                    # just the sensitivity evaluation

        else:  # this is if it's OpenMDAO
            for i, iter_type in enumerate(self.iter_type):
                key = f"{self.solver_name}|{i + 1}"  # OpenMDAO uses 1-indexing
                if i in self.deriv_keys:
                    self.iter_type[i] = 1.0

            # If no derivative info is saved, we don't know which iterations are major.
            # Treat all iterations as major.
            if len(self.deriv_keys) < 1:
                self.iter_type[:] = 1.0

    def SaveDBData(self, db, data_all, data_major, OpenMDAO, data_str):
        """Method to save the information within the database corresponding
        to a certain key to the relevant dictionaries within the Display
        object. This method is called twice, once for the design variables
        and the other for the outputs."""

        # Loop over each optimization iteration
        for i, iter_type in enumerate(self.iter_type):
            # If this is an OpenMDAO file, the keys are of the format
            # 'rank0:SNOPT|1', etc
            if OpenMDAO:
                key = f"{self.solver_name}|{i + 1}"  # OpenMDAO uses 1-indexing
            else:  # Otherwise the keys are simply a number
                key = "%d" % i

            # Do this for both major and minor iterations
            if self.iter_type[i]:
                # Get just the info in the dict for this iteration
                iter_data = db[key][data_str]

                # Loop through each key within this iteration
                for key in sorted(iter_data):
                    # Format a new_key string where we append a modifier
                    # if we have multiple history files
                    new_key = key + f"{self.histIndex}"

                    # If this key is not in the data dictionaries, add it
                    if new_key not in data_all:
                        data_all[new_key] = []
                        data_major[new_key] = []

                    # Process the data from the key. Convert it to a numpy
                    # array, keep only the real part, squeeze any 1-dim
                    # axes out of it, then flatten it.
                    data = np.squeeze(np.array(iter_data[key]).real).flatten()

                    # Append the data to the entries within the dictionaries.
                    data_all[new_key].append(data)
                    if self.iter_type[i] == 1:
                        data_major[new_key].append(data)

    def SaveOpenMDAOData(self, db):
        """Examine the OpenMDAO dict and save tags if the variables are
        objectives (o), constraints (c), or design variables (dv)."""

        # Loop over each key in the metadata db
        for tag in db:
            # Only look at variables and unknowns
            if tag in ["Unknowns", "Parameters"]:
                for old_item in db[tag]:
                    # We'll rename each item, so we need to get the old item
                    # name and modify it
                    item = old_item + f"{self.histIndex}"

                    # Here we just have an open parenthesis, and then we will
                    # add o, c, or dv. Note that we could add multiple flags
                    # to a single item. That's why we have a sort of convoluted
                    # process of adding the tags.
                    new_key = item + " ("
                    flag_list = []

                    # Check each flag and see if they have the relevant entries
                    # within the dict; if so, tag them.
                    for flag in db[tag][old_item]:
                        if "is_objective" in flag:
                            flag_list.append("o")
                        if "is_desvar" in flag:
                            flag_list.append("dv")
                        if "is_constraint" in flag:
                            flag_list.append("c")

                    # Create the new_key based on the flags for each variable
                    for flag in flag_list:
                        if flag == flag_list[-1]:
                            new_key += flag + ")"
                        else:
                            new_key += flag + ", "

                    # If there are actually flags to add, pop out the old items
                    # in the dict and re-add them with the new name.
                    if flag_list:
                        try:
                            if "dv" in flag_list:
                                self.var_data_all[new_key] = self.func_data_all.pop(item)
                                self.var_data_major[new_key] = self.func_data_major.pop(item)

                            else:
                                self.func_data_all[new_key] = self.func_data_all.pop(item)
                                self.func_data_major[new_key] = self.func_data_major.pop(item)

                        except KeyError:
                            pass
