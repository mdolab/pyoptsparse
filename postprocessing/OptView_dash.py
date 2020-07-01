#!/usr/bin/python
import dash
import dash_core_components as dcc
import dash_html_components as html
import plotly.graph_objs as go
from plotly import tools
import numpy as np
import argparse
from sqlitedict import SqliteDict
import shelve
import sys

from OptView_baseclass import OVBaseClass


major_python_version = sys.version_info[0]
parser = argparse.ArgumentParser()
parser.add_argument(
    "histFile", nargs="*", type=str, default="opt_hist.hst", help="Specify the history file to be plotted"
)

args = parser.parse_args()
histList = args.histFile


class ReadOptHist(OVBaseClass):
    """
    Container for display parameters, properties, and objects.
    This includes a canvas for MPL plots and a bottom area with widgets.
    """

    def __init__(self, histList):

        # Initialize lists, dicts, and save inputs from user
        self.arr_active = 0
        self.plots = []
        self.annotate = None
        self.histList = histList
        self.bounds = {}
        self.scaling = {}
        self.color_bounds = [0.0, 0.0]

        # Actually setup and run the GUI
        self.OptimizationHistory()

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
            except:  # Bare except because error is not in standard Python. # noqa: E722
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
                if major_python_version == 3:
                    for string in db.keys():
                        string = string.split("|")
                else:
                    string = db.keys()[-1].split("|")

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
                self.iter_type = np.zeros(nkey)

                # Check to see if there is bounds information in the db file.
                # If so, add them to self.bounds to plot later.
                try:
                    info_dict = db["varInfo"].copy()
                    info_dict.update(db["conInfo"])
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
                        scaling_dict[key + histIndex] = info_dict[key]["scale"]
                    self.bounds.update(bounds_dict)
                    self.scaling.update(scaling_dict)
                except KeyError:
                    pass

                # Check to see if there is proper saved info about iter type
                if "iu0" in db["0"].keys():
                    if db[db["last"]]["iu0"] > 0:
                        self.storedIters = True
                    else:
                        self.storedIters = False
                else:
                    self.storedIters = False

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
            # Loop over each optimization iteration
            for i, iter_type in enumerate(self.iter_type):

                # If this is an OpenMDAO file, the keys are of the format
                # 'rank0:SNOPT|1', etc
                key = "%d" % i

                # Only actual optimization iterations have 'funcs' in them.
                # pyOptSparse saves info for two iterations for every
                # actual major iteration. In particular, one has funcs
                # and the next has funcsSens, but they're both part of the
                # same major iteration.
                if any("funcs" == s for s in db[key].keys()):

                    # If this iteration has 'funcs' within it, but it's not
                    # flagged as major, then it's a minor iteration.
                    if i == 0:
                        self.iter_type[i] = 1
                    else:
                        self.iter_type[i] = 2

                    if not self.storedIters:  # Otherwise, use a spotty heuristic to see if the
                        # iteration is major or not. NOTE: this is often
                        # inaccurate, especially if the optimization used
                        # gradient-enhanced line searches.
                        try:
                            keyp1 = "%d" % (i + 1)
                            db[keyp1]["funcsSens"]
                            self.iter_type[i] = 1  # for 'major' iterations
                        except KeyError:
                            self.iter_type[i] = 2  # for 'minor' iterations

                else:
                    self.iter_type[i] = 0  # this is not a real iteration,
                    # just the sensitivity evaluation

            # Loop over each optimization iteration
            for i, iter_type in enumerate(self.iter_type):

                if iter_type == 0:
                    continue

                key = "%d" % i

        else:  # this is if it's OpenMDAO
            for i, iter_type in enumerate(self.iter_type):
                key = "{}|{}".format(self.solver_name, i + 1)  # OpenMDAO uses 1-indexing
                if i in self.deriv_keys:
                    self.iter_type[i] = 1.0

            # If no derivative info is saved, we don't know which iterations are major.
            # Treat all iterations as major.
            if len(self.deriv_keys) < 1:
                self.iter_type[:] = 1.0

    def SaveDBData(self, db, data_all, data_major, OpenMDAO, data_str):
        """ Method to save the information within the database corresponding
            to a certain key to the relevant dictionaries within the Display
            object. This method is called twice, once for the design variables
            and the other for the outputs. """

        # Loop over each optimization iteration
        for i, iter_type in enumerate(self.iter_type):

            # If this is an OpenMDAO file, the keys are of the format
            # 'rank0:SNOPT|1', etc
            if OpenMDAO:
                key = "{}|{}".format(self.solver_name, i + 1)  # OpenMDAO uses 1-indexing
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
                    new_key = key + "{}".format(self.histIndex)

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
        """ Examine the OpenMDAO dict and save tags if the variables are
            objectives (o), constraints (c), or design variables (dv). """

        # Loop over each key in the metadata db
        for tag in db:

            # Only look at variables and unknowns
            if tag in ["Unknowns", "Parameters"]:
                for old_item in db[tag]:

                    # We'll rename each item, so we need to get the old item
                    # name and modify it
                    item = old_item + "{}".format(self.histIndex)

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

    def refresh_history(self):
        """
        Refresh opt_his data if the history file has been updated.
        """
        # old_funcs = []
        # for key in self.func_data:
        #     old_funcs.append(key)
        # old_vars = []
        # for key in self.var_data:
        #     old_vars.append(key)

        self.OptimizationHistory()

        # new_funcs = []
        # for key in self.func_data:
        #     new_funcs.append(key)
        # new_vars = []
        # for key in self.var_data:
        #     new_vars.append(key)
        #
        # if not (old_funcs == new_funcs and old_vars == new_vars):
        #     self.var_search('dummy')


Opt = ReadOptHist(histList)
app = dash.Dash(__name__)
app.layout = html.Div(
    [
        html.Div(
            className="banner",
            children=[
                # Change App Name here
                html.Div(
                    className="container scalable",
                    children=[
                        # Change App Name here
                        html.H2(html.A("MACH Dashboard", style={"text-decoration": "none", "color": "inherit"})),
                    ],
                ),
            ],
        ),
        html.Div(
            id="body",
            className="container scalable",
            children=[
                html.Div(
                    className="row",
                    children=[
                        html.Div(
                            className="ten columns",
                            style={
                                "min-width": "24.5%",
                                "height": "calc(100vh - 120px)",
                                "margin-top": "5px",
                                # Remove possibility to select the text for better UX
                                "user-select": "none",
                                "-moz-user-select": "none",
                                "-webkit-user-select": "none",
                                "-ms-user-select": "none",
                            },
                            children=[
                                html.Div(
                                    [
                                        dcc.Dropdown(
                                            id="dvar",
                                            options=[{"label": i, "value": i} for i in Opt.var_data_all.keys()],
                                            multi=True,
                                            placeholder="Choose design group",
                                        ),
                                    ],
                                    style={"width": "50%", "display": "inline-block"},
                                ),
                                html.Div(
                                    [
                                        dcc.Dropdown(
                                            id="func",
                                            options=[{"label": i, "value": i} for i in Opt.func_data_all.keys()],
                                            multi=True,
                                            placeholder="Choose func group",
                                        ),
                                    ],
                                    style={"width": "50%", "display": "inline-block"},
                                ),
                                html.Div(
                                    [dcc.Dropdown(id="dvar-child", multi=True, placeholder="Choose design var")],
                                    style={"width": "50%", "display": "inline-block"},
                                ),
                                html.Div(
                                    dcc.Dropdown(id="func-child", multi=True, placeholder="Choose func var"),
                                    style={"width": "50%", "display": "inline-block"},
                                ),
                                html.Div(
                                    dcc.RadioItems(
                                        id="axis_scale",
                                        options=[
                                            {"label": "linear", "value": "linear"},
                                            {"label": "log", "value": "log"},
                                        ],
                                        value="linear",
                                        labelStyle={"display": "inline-block"},
                                    ),
                                    style={"width": "50%", "display": "inline-block"},
                                ),
                                html.Div(
                                    dcc.RadioItems(
                                        id="scale_type",
                                        options=[
                                            {"label": "shared axes", "value": "shared"},
                                            {"label": "multiple axes", "value": "multi"},
                                        ],
                                        value="shared",
                                        labelStyle={"display": "inline-block"},
                                    ),
                                    style={"width": "50%", "display": "inline-block"},
                                ),
                                dcc.Graph(id="2Dscatter"),
                                dcc.Interval(
                                    id="interval-component", interval=1 * 1000, n_intervals=0  # in milliseconds
                                ),
                                html.Div(id="hidden-div", style={"display": "none"}),
                            ],
                        )
                    ],
                )
            ],
        ),
    ]
)

external_css = [
    # Normalize the CSS
    "https://cdnjs.cloudflare.com/ajax/libs/normalize/7.0.0/normalize.min.css",
    # Fonts
    "https://fonts.googleapis.com/css?family=Open+Sans|Roboto",
    "https://maxcdn.bootstrapcdn.com/font-awesome/4.7.0/css/font-awesome.min.css",
]

for css in external_css:
    app.css.append_css({"external_url": css})

app.config.suppress_callback_exceptions = True


@app.callback(
    dash.dependencies.Output("2Dscatter", "figure"),
    [
        dash.dependencies.Input("dvar-child", "value"),
        dash.dependencies.Input("func-child", "value"),
        dash.dependencies.Input("axis_scale", "value"),
        dash.dependencies.Input("scale_type", "value"),
        dash.dependencies.Input("hidden-div", "value"),
    ],
)
def update_2d_scatter(dvar, func, axisscale, type, n):
    trace = []
    i = 0
    if dvar:
        for var in dvar:
            index = var.split("_")[-1]
            varname = var[::-1].replace(index + "_", "", 1)[::-1]
            trace.append(
                go.Scatter(
                    x=range(Opt.num_iter),
                    y=[data[int(index)] for data in Opt.var_data[varname]],
                    name=var,
                    mode="lines+markers",
                )
            )
            i += 1

    if func:
        for var in func:
            index = var.split("_")[-1]
            varname = var[::-1].replace(index + "_", "", 1)[::-1]
            trace.append(
                go.Scatter(
                    x=range(Opt.num_iter),
                    y=[data[int(index)] for data in Opt.func_data_all[varname]],
                    name=var,
                    mode="lines+markers",
                )
            )
            i += 1

    fig = {}
    fig["layout"] = {}
    if dvar or func:
        if type == "multi":
            fig = tools.make_subplots(rows=i, cols=1)
            for k in range(i):
                fig.append_trace(trace[k], k + 1, 1)
        else:
            fig["data"] = trace

    fig["layout"].update(
        xaxis={
            "title": {
                "text": "iterations",
                "font": {"family": "Courier New, monospace", "size": 24, "color": "#7f7f7f"},
            },
            "type": axisscale,
        },
        yaxis={
            "title": {"text": "Data", "font": {"family": "Courier New, monospace", "size": 24, "color": "#7f7f7f"}},
            "type": axisscale,
        },
        height=900,
        showlegend=True,
    )

    return fig


@app.callback(dash.dependencies.Output("dvar-child", "options"), [dash.dependencies.Input("dvar", "value")])
def update_dvar_child(dvar):
    strlist = []
    if dvar:
        for var in dvar:
            num = len(Opt.var_data[var][0])
            strlist += [var + "_" + str(i) for i in range(num)]
    return [{"label": i, "value": i} for i in strlist]


@app.callback(dash.dependencies.Output("func-child", "options"), [dash.dependencies.Input("func", "value")])
def update_func_child(func):
    strlist = []
    if func:
        for var in func:
            num = len(Opt.func_data_all[var][0])
            strlist += [var + "_" + str(i) for i in range(num)]

    return [{"label": i, "value": i} for i in strlist]


@app.callback(
    dash.dependencies.Output("hidden-div", "value"), [dash.dependencies.Input("interval-component", "n_intervals")]
)
def update_opt_history(n):
    Opt.refresh_history()
    return n


if __name__ == "__main__":
    app.run_server(debug=True, port=50002, host="0.0.0.0")
