"""
This dash version of OptView makes use of the new API
to read in the history file, rather than using the
OptView baseclass. This should be more maintainable
for adding new features or displaying new information with OptView.
"""

# Standard Python modules
import argparse
import json

# External modules
import dash
import dash_core_components as dcc
import dash_html_components as html
import numpy as np
from plotly import graph_objs as go
from plotly import subplots

# First party modules
from pyoptsparse import History

# Read in the history files given by user
parser = argparse.ArgumentParser()
parser.add_argument(
    "histFile", nargs="*", type=str, default="opt_hist.hst", help="Specify the history file to be plotted"
)
args = parser.parse_args()

# List of strings (the history file names)
histListArgs = args.histFile

# Create list of associated labels (A, B, C, etc) for each history file if there are more than one
index = ord("A")
fileLabels = ["" for file in histListArgs]
if len(fileLabels) > 1:
    for i in range(len(fileLabels)):
        fileLabels[i] = chr(index)
        index = (index + 1 - 65) % 26 + 65

# Color scheme for graphing traces
colors = ["#636EFA", "#EF553B", "#00CC96", "#AB63FA", "#FFA15A", "#19D3F3", "#FF6692", "#B6E880", "#FF97FF", "#FECB52"]


# ====================================================================================
#  Defining dash app & layout
# ====================================================================================

app = dash.Dash(__name__)
# Override exceptions for when elements are defined without initial input
app.config.suppress_callback_exceptions = True

app.layout = html.Div(
    children=[
        html.Nav(
            children=[
                html.Img(
                    src=app.get_asset_url("/logo.png"), style={"height": "4rem", "width": "13rem", "margin": "0.2rem"}
                ),
                html.Div(["OptView"], style={"fontWeight": "550", "fontSize": "3.0rem", "color": "#193D72"}),
            ],
            style={"display": "flex", "backgroundColor": "#EEEEEE"},
        ),
        html.Div(
            [
                html.Div(
                    [
                        html.Div(
                            [
                                html.Div(
                                    [
                                        html.Div(
                                            ["History Files"],
                                            style={"color": "#193D72", "fontWeight": "normal", "fontSize": "15px"},
                                        ),
                                        html.Ul(
                                            [
                                                html.Li(
                                                    [fileLabels[i] + " " + histListArgs[i]], style={"fontSize": "12px"}
                                                )
                                                for i in range(len(histListArgs))
                                            ]
                                        ),
                                    ],
                                ),
                                html.Div(
                                    [
                                        html.Div(
                                            ["Design Groups"],
                                            className="dvarGroup",
                                            style={"color": "#193D72", "fontWeight": "normal", "fontSize": "1.5rem"},
                                        ),
                                        dcc.Dropdown(
                                            id="dvarGroup",
                                            placeholder="Select design group(s)...",
                                            multi=True,
                                            style={"color": "#676464", "fontSize": "1.2rem"},
                                        ),
                                    ],
                                ),
                                html.Div(
                                    [
                                        html.Div(
                                            ["Design Variables"],
                                            style={"color": "#193D72", "fontWeight": "normal", "fontSize": "1.5rem"},
                                        ),
                                        dcc.Dropdown(
                                            id="dvarChild",
                                            placeholder="Select design var(s)...",
                                            multi=True,
                                            style={"color": "#676464", "fontSize": "1.2rem"},
                                        ),
                                    ],
                                ),
                            ],
                            style={"marginRight": "1.0rem", "width": "20.0rem"},
                        ),
                        html.Div(
                            [
                                html.Div(
                                    [
                                        html.Div(
                                            ["Function Groups"],
                                            style={"color": "#193D72", "fontWeight": "normal", "fontSize": "1.5rem"},
                                        ),
                                        dcc.Dropdown(
                                            id="funcGroup",
                                            placeholder="Select function group(s)...",
                                            multi=True,
                                            style={"color": "#676464", "fontSize": "1.2rem"},
                                        ),
                                    ],
                                ),
                                html.Div(
                                    [
                                        html.Div(
                                            ["Function Variables"],
                                            style={"color": "#193D72", "fontWeight": "normal", "fontSize": "1.5rem"},
                                        ),
                                        dcc.Dropdown(
                                            id="funcChild",
                                            placeholder="Select function var(s)...",
                                            multi=True,
                                            style={"color": "#676464", "fontSize": "1.2rem"},
                                        ),
                                    ],
                                ),
                            ],
                            style={"marginRight": "3.0rem", "width": "20.0rem"},
                        ),
                        html.Div(
                            [
                                html.Div(
                                    [
                                        html.Div(
                                            ["Optimization Groups"],
                                            style={"color": "#193D72", "fontWeight": "normal", "fontSize": "1.5rem"},
                                        ),
                                        dcc.Dropdown(
                                            id="optGroup",
                                            placeholder="Select optimization group(s)...",
                                            multi=True,
                                            style={"color": "#676464", "fontSize": "1.2rem"},
                                        ),
                                    ],
                                ),
                                html.Div(
                                    [
                                        html.Div(
                                            ["Optimization Variables"],
                                            style={"color": "#193D72", "fontWeight": "normal", "fontSize": "1.5rem"},
                                        ),
                                        dcc.Dropdown(
                                            id="optChild",
                                            placeholder="Select optimization var(s)...",
                                            multi=True,
                                            style={"color": "#676464", "fontSize": "1.2rem"},
                                        ),
                                    ],
                                ),
                            ],
                            style={"marginRight": "3.0rem", "width": "20.0rem"},
                        ),
                        html.Div(
                            [
                                html.Div(
                                    ["Graph"], style={"color": "#193D72", "fontWeight": "normal", "fontSize": "1.5rem"}
                                ),
                                dcc.RadioItems(
                                    id="plot-type",
                                    options=[
                                        {"label": "Stacked", "value": "stacked"},
                                        {"label": "Shared", "value": "shared"},
                                    ],
                                    value="shared",
                                    style={"color": "#676464", "fontSize": "1.2rem"},
                                ),
                            ],
                            style={"width": "30x", "marginRight": "2%"},
                        ),
                        html.Div(
                            [
                                html.Div(
                                    ["Additional"],
                                    style={"color": "#193D72", "fontWeight": "normal", "fontSize": "1.5rem"},
                                ),
                                dcc.Checklist(
                                    id="data-type",
                                    options=[
                                        {"label": "Show Bounds", "value": "bounds"},
                                        {"label": "Show Major Iterations", "value": "major"},
                                        {"label": "Log Scale", "value": "log"},
                                        {"label": "Apply Scaling Factor", "value": "scale"},
                                        {"label": "Show Min/Max of Group(s)", "value": "minMax"},
                                        {"label": "Show Absolute Delta", "value": "delta"},
                                        {"label": "Auto-Refresh (Toggle on to enter Rate)", "value": "refresh"},
                                    ],
                                    style={"color": "#676464", "fontSize": "1.2rem"},
                                    value=["major"],
                                ),
                                dcc.Input(
                                    id="refresh-rate",
                                    type="number",
                                    placeholder="Refresh Rate (s)",
                                    style={"fontSize": "1.2rem", "marginLeft": "2rem"},
                                ),
                                # dcc.Input(
                                #     id='figure-export-width',
                                #     type='number',
                                #     placeholder='Enter PX Width',
                                #     style={"fontSize":"1.2rem"}
                                # ),
                                # dcc.Input(
                                #     id='figure-export-height',
                                #     type='number',
                                #     placeholder='Enter PX Height',
                                #     style={"fontSize":"1.2rem"}
                                # )
                            ],
                        ),
                    ],
                    style={"display": "flex", "flexDirection": "column", "padding": "0.5rem"},
                ),
                html.Div(
                    [
                        dcc.Graph(
                            id="plot",
                            config={
                                "scrollZoom": True,
                                "showTips": True,
                                "toImageButtonOptions": {
                                    "format": "png",
                                    "filename": "custom_image",
                                    "height": 500,
                                    "width": 700,
                                    "scale": 10,
                                },
                            },
                        )
                    ],
                    style={"width": "100%"},
                ),
                dcc.Interval(
                    id="interval-component",
                    interval=1 * 1000,
                    n_intervals=0,  # in milliseconds
                    # max_intervals=5
                ),
                html.Div(id="hidden-div", style={"display": "none"}),
            ],
            style={"display": "flex"},
        ),
    ],
)

# ====================================================================================
#  Helper functions for Dash Callback Functions found further down
# ====================================================================================


def getValues(name, dataType, hist):
    """
    Uses user input and returns a data dictionary with values for the requested value of interest, by calling the getValues() method
    from pyOpt_history.py, adjusting the requested values using the dataType specifications from the user.

    Parameters
    ----------
    name : str
        The value of interest, can be the name of any DV, objective or constraint that a user selects

    dataType : list of str
        Contains dataType str values selected by user (i.e scale, major, delta)

    hist: History object
        This is created with the history file for the given value of interest using pyOpt-history.py

    Returns
    -------
    data dictionary
        Contains the iteration data for the value requested. It returns a data dictionary with the key
        as the 'name' requested and the value a numpy array where with the first dimension equal to the
        number of iterations requested (depending on if major is in dataType or not).

    Example
    _______
    getValues(name="xvars", dataType=[['major']], hist):
        where the xvars group contains 3 variables and 5 major iterations, would return an array with the first dimension
        being 5, the number of iterations, and the second dimension being 3, where there are 3 values for the 3 variables per iteration.
    """

    if dataType:
        values = hist.getValues(names=name, major=("major" in dataType), scale=("scale" in dataType))
        if "delta" in dataType:
            tempValues = values[name].copy()
            for i in list(range(len(values[name]))):
                for j in list(range(len(values[name][i]))):
                    if i != 0:
                        values[name][i][j] = abs(values[name][i][j] - tempValues[i - 1][j])
                    else:
                        values[name][i][j] = 0

    else:
        values = hist.getValues(names=name, major=False)
    return values


def addVarTraces(var, trace, dataType, names, groupType):
    """
    Adds traces to the plotly figure trace[] list for the passed in variable, based on specified user input including data tpye
    and plot type. Also adds upper and lower bounds traces if specified by user in dataType.

    If 'bounds' is in dataType, three traces are added: the var's trace, and the upper and lower bound traces for the var.
    If not, one trace is added, being the passed in var's trace.

    Parameters
    ----------
    var : str
        The name of the variable to add a trace for
        If multiple history files are being used - Formatted as: OPT/DV/FUNC-GROUPNAME_HIST-FILE-LABEL-INDEX_VAR-INDEX
            EX: xvars_A_0 -> This would be the first variable in the xvars group from file A
        One history file:
            EX: EX: xvars_0 -> This would be the first variable in the xvars group
        *NOTE* If the xvars group only has one scalar quantity, the corresponding name will only be xvars with no _INDEX

    trace: list
        List containing the traces to be displayed on the plotly figure.

    dataType : list of str
        Contains dataType str values selected by user (i.e scale, major, delta, bounds)

    names: list of str
        List of either the function, optimization, or   DV group names that have been selected by the user.

    groupType: str
        'dvar', 'func', or 'opt' depending which group the specified var is in

    Returns
    -------
    Nothing
    """
    # A list of the History objects for each history file argument
    histList = [History(fileName) for fileName in histListArgs]
    var = str(var)
    # Initializing History object associated with the passed in var
    hist = histList[0]
    # Variable index (EX: var=xvars_0 -> indexVar=0)
    indexVar = 0
    # History file index (EX: var=xvars_A_0 -> indexHist=A)
    indexHist = "A"

    # # Group name for the passed in variable (EX: var=xvars_A_0 -> varGroupName=xvars_A)
    # varGroupName = var
    # Variable key name that is used in history file (removed _INDEX or _indexHist) (EX: var=xvars_0/xvars_A_0  -> varName=xvars)
    varName = var

    # Set varName, indexHist, and varGroup name accordingly based on # of history files and # of variables in var's group
    if var in names:
        # If the var is in the names list, then it does not have an _INDEX appointed and must be apart of a group with one qty
        if len(histList) > 1:
            indexHist = var.split("_")[-1]
            varName = var[::-1].replace(indexHist + "_", "", 1)[::-1]
    else:
        indexVar = var.split("_")[-1]
        reversedIndexVar = indexVar[::-1]
        # varGroupName = var[::-1].replace(indexVar + "_", "", 1)[::-1]
        if len(histList) > 1:
            indexHist = var.split("_")[-2]
            varName = var[::-1].replace(reversedIndexVar + "_" + indexHist + "_", "", 1)[::-1]
        else:
            varName = var[::-1].replace(reversedIndexVar + "_", "", 1)[::-1]

    # Set hist as the appropriate History object for the var
    hist = histList[ord(indexHist) % 65]
    # Dictionary with DV or Func group info for the var, with keys 'scale', 'lower', 'upper',
    # where each key corresponds to a numpy array of size # variables, with each index's value corresponding to that indexed variable in the group
    info = hist.getDVInfo() if groupType == "dvar" else hist.getConInfo()
    # Needed trace data based on var and dataType needed
    data = getValues(name=varName, dataType=dataType, hist=hist)

    # Add trace for var to trace[] list
    trace.append(
        go.Scatter(
            x=list(range(len(data[varName]))),
            y=[data.real[int(indexVar)] for data in data[varName]],
            name=var,
            marker_color=colors[(len(trace) - 1) % len(colors)],
            mode="lines+markers",
            marker={"size": 3},
            line={"width": 1},
        )
    )
    # Add upper & lower bounds traces if 'bounds' was selected
    if dataType and "bounds" in dataType and varName != "obj" and groupType != "opt":
        # Create lower + upper bounds array to plot, and scale if scaling is needed. If no 'lower' or 'upper' information, leave trace blank
        if "scale" in dataType:
            # HACK/TODO: Used np.print_1d to convert all ['scale'] values to a 1d array ('scale' gives a scalar qty for groups with only one variable)
            scaleData = np.atleast_1d(info[varName]["scale"])
            # HACK/TODO: ['scale'] array should be the size of the variables in a group but some seem to
            # be one constant, so this checks if it should be made as an array of the size of the vars
            if var not in names and len(scaleData) == 1:
                scaleData = [info[varName]["scale"]] * len(data[varName])
            scaleFactor = np.atleast_1d(info[varName]["scale"])[int(indexVar)].real
            lowerB = (
                [info[varName]["lower"][int(indexVar)] * scaleFactor] * len(data[varName])
                if info[varName]["lower"][0]
                else []
            )
            upperB = (
                [info[varName]["upper"][int(indexVar)] * scaleFactor] * len(data[varName])
                if info[varName]["upper"][0]
                else []
            )
        else:
            lowerB = [info[varName]["lower"][int(indexVar)]] * len(data[varName])
            upperB = [info[varName]["upper"][int(indexVar)]] * len(data[varName])
        # Add lower + upper bound traces to trace list
        trace.append(
            go.Scatter(
                x=list(range(len(data[varName]))),
                y=lowerB,
                name=var + "_LB",
                marker_color=colors[(len(trace) - 2) % len(colors)],
                mode="lines",
                line={"dash": "dash", "width": 2},
            )
        )
        trace.append(
            go.Scatter(
                x=list(range(len(data[varName]))),
                y=upperB,
                name=var + "_UB",
                marker_color=colors[(len(trace) - 3) % len(colors)],
                mode="lines",
                line={"dash": "dash", "width": 2},
            )
        )


def addMaxMinTraces(var, trace, dataType, histList):
    """
    Adds min and max to the plotly figure trace[] list for the passed in variable group.

    Parameters
    ----------
    var : str
        The name of the variable group to add min/max traces for.
            EX: xvars

    trace: list
        List containing the traces to be displayed on the plotly figure.

    dataType : list of str
        Contains dataType str values selected by user (i.e scale, major, delta, bounds)

    histList: list of str
        List of str file names to be turned into history objects.

    Returns
    -------
    Nothing
    """
    var = str(var)
    # Variable key name that is used in history file (removed _INDEX or _indexHist) (EX: var=xvars_0/xvars_A_0  -> varName=xvars)
    varName = var
    hist = histList[0]
    if len(histList) > 1:
        indexHist = var.split("_")[-1]
        hist = histList[ord(indexHist) % 65]
        varName = var[::-1].replace(indexHist + "_", "", 1)[::-1]
    data = getValues(name=varName, dataType=dataType, hist=hist)
    # Append minMax traces
    # if(dataType and 'minMax' in dataType):
    trace.append(
        go.Scatter(
            x=list(range(len(data[varName]))),
            y=[max(arr.real) for arr in data[varName]],
            name=var + "_max",
            mode="lines+markers",
            marker={"size": 3},
            line={"width": 1},
        )
    )
    trace.append(
        go.Scatter(
            x=list(range(len(data[varName]))),
            y=[min(arr.real) for arr in data[varName]],
            name=var + "_min",
            mode="lines+markers",
            marker={"size": 3},
            line={"width": 1},
        )
    )


# ====================================================================================
#  Dash Callback Functions (Triggered by front-end interface input changes)
# ====================================================================================


@app.callback(
    dash.dependencies.Output("interval-component", "max_intervals"), [dash.dependencies.Input("data-type", "value")]
)
def update_refresh_intervals(dataType):
    """
    This determines if the interval component should run infinitely based on if the autorefresh
    setting is turned on.

    Parameters
    ----------
    dataType : list of str
        Contains dataType str values selected by user (i.e scale, major, delta, bounds)

    Returns
    -------
    int
        -1 if there should be an infinite amount of intervals, 0 if it should not autorefresh
    """
    if dataType and "refresh" in dataType:
        return -1
    else:
        return 0


@app.callback(
    dash.dependencies.Output("interval-component", "interval"), [dash.dependencies.Input("refresh-rate", "value")]
)
def update_refresh_rate(rate):
    """
    This determines the rate that interval component  should run at in ms based on user input.

    Parameters
    ----------
    rate: int
        The rate inputted by the user in seconds.

    Returns
    -------
    int
        The rate in ms for the interval component to run at.
    """
    if rate:
        return rate * 1000
    else:
        return 10000


@app.callback(dash.dependencies.Output("refresh-rate", "disabled"), [dash.dependencies.Input("data-type", "value")])
def update_refresh_input(dataType):
    """
    This prevents the user from inputting the refresh rate if the auto refresh feature is not turned on.

    Parameters
    ----------
     dataType : list of str
        Contains dataType str values selected by user (i.e scale, major, delta, bounds)

    Returns
    -------
    bool
        True to disable the input box, False to disable it
    """
    if dataType and "refresh" in dataType:
        return False
    else:
        return True


@app.callback(
    dash.dependencies.Output("hidden-div", "children"),
    [
        dash.dependencies.Input("interval-component", "n_intervals"),
        dash.dependencies.Input("interval-component", "max_intervals"),
    ],
)
def update_opt_history(n, max):
    """
    This saves the updated function and design group names in an object in the hidden-div, that updates
    based on the interval-component's number of intervals that runs on interval

    Parameters
    ----------
     n: int
        The current number of intervals ran.
        Not explicitly used in this function, but kept as an input dependency so this function is called on interval.

    max: int
        The max number of intervals the interval-component runs on.
        Not explicitly used in this function, but kept as an input dependency so that if auto refresh is not on,
        this function is called once when the max_intervals is set/updated as 0.

    Returns
    -------
    dictionary
        Saves a dictionary of the 'dvNames' and 'funcNames' for the groups as the history file refreshes. If autorefresh is not
        on, this will be called once.

        If multiple history files are being used - Formatted as: DV/FUNC-GROUPNAME_HIST-FILE-LABEL-INDEX
          EX: xvars_A -> This would be the the design group 'xvars' from the file indexed A
    """
    # Re-reads in history list data and passes to History API
    try:
        # If AutoRefresh is on and the data has not yet loaded, this line will throw an error
        histList = [History(fileName) for fileName in histListArgs]
    except Exception:
        print("History file data is not processed yet")
        return {}
    else:
        # Re-Saving names for drop down menus
        index = ord("A")
        dvNames = []
        funcNames = []
        optNames = []
        if len(histList) > 1:
            for hist in histList:
                dvNames += [name + "_" + chr(index) for name in hist.getDVNames()]
                conNames = [name + "_" + chr(index) for name in hist.getConNames()]
                objNames = [name + "_" + chr(index) for name in hist.getObjNames()]
                funcNames += conNames + objNames
                optNames = list(
                    set(histList[0].getIterKeys()).difference({"xuser", "funcs", "fail", "nMajor", "nMinor", "isMajor"})
                )
                optNames = [name + "_" + chr(index) for name in optNames]
                index = (index + 1 - 65) % 26 + 65
        else:
            hist = histList[0]
            dvNames = hist.getDVNames()
            conNames = hist.getConNames()
            objNames = hist.getObjNames()
            funcNames = conNames + objNames
            optNames = list(
                set(histList[0].getIterKeys()).difference({"xuser", "funcs", "fail", "nMajor", "nMinor", "isMajor"})
            )

        historyInfo = {
            "dvNames": dvNames,
            "funcNames": funcNames,
            "optNames": optNames
            # 'values' : [hist.getValues(major=False, scale=False) for hist in histList],
            # 'valuesMajor' : [hist.getValues(major=True) for hist in histList],
            # 'valuesScale' : [hist.getValues(scale=True, major=False) for hist in histList],
            # 'valuesScaleMajor' : [hist.getValues(scale=True, major=True) for hist in histList],
        }

        # Parses object into JSON format to return to the hidden-div.
        # Returns all needed history object info to div
        # obj info, con info, dv info, getValues with: major=true, major=false, major=true+scale=true, major=false+scale=true
        # print(len(historyInfo['values']))
        return json.dumps(historyInfo)


@app.callback(dash.dependencies.Output("funcGroup", "options"), [dash.dependencies.Input("hidden-div", "children")])
def update_funcGroupOptions(historyInfo):
    """
    This populates the 'Function Group' dropdown based on the info from the dictionary saved in the
    'hidden-div' component.

    Parameters
    ----------
    historyInfo: dictionary
        This is a dictionary with the keys 'dvNames', 'funcNames', and 'optNames' holding a list each of the group names for the Group dropdowns.

    Returns
    -------
    list of dictionaries
        List of dictionaries containing the options for the 'Function Group' dropdown.

        If multiple history files are being used - Formatted as: OPT/DV/FUNC-GROUPNAME_HIST-FILE-LABEL-INDEX
          EX: twist_B -> This would be the the design group 'xvars' from the file indexed B
    """
    if historyInfo:
        # If historyInfo isn't null, then it is ready to be parsed to update the function group dropdown options
        historyInfo = json.loads(historyInfo)
        options = [{"label": i, "value": i} for i in historyInfo["funcNames"]]
        return options
    else:
        # If historyInfo is null, keep the options as empty
        return []


@app.callback(dash.dependencies.Output("dvarGroup", "options"), [dash.dependencies.Input("hidden-div", "children")])
def update_dvarGroupOptions(historyInfo):
    """
    This populates the 'Design Group' dropdown based on the info from the dictionary saved in the
    'hidden-div' component.

    Parameters
    ----------
    historyInfo: dictionary
        This is a dictionary with the keys 'dvNames', 'funcNames', and 'optNames' holding a list each of the group names for the Group dropdowns.

    Returns
    -------
    list of dictionaries
        List of dictionaries containing the options for the 'Design Group' dropdown.

    If multiple history files are being used - Formatted as: OPT/DV/FUNC-GROUPNAME_HIST-FILE-LABEL-INDEX
          EX: twist_B -> This would be the the design group 'xvars' from the file indexed B
    """
    if historyInfo:
        # If historyInfo isn't null, then it is ready to be parsed to update the function group dropdown options
        historyInfo = json.loads(historyInfo)
        options = [{"label": i, "value": i} for i in historyInfo["dvNames"]]
        return options
    else:
        # If historyInfo is null, keep the options as empty
        return []


@app.callback(dash.dependencies.Output("optGroup", "options"), [dash.dependencies.Input("hidden-div", "children")])
def update_optGroupOptions(historyInfo):
    """
    This populates the 'Optimization Group' dropdown based on the info from the dictionary saved in the
    'hidden-div' component.

    Parameters
    ----------
    historyInfo: dictionary
        This is a dictionary with the keys 'dvNames', 'funcNames', and 'optNames' holding a list each of the group names for the Group dropdowns.

    Returns
    -------
    list of dictionaries
        List of dictionaries containing the options for the 'Optimization Group' dropdown.

    If multiple history files are being used - Formatted as: OPT/DV/FUNC-GROUPNAME_HIST-FILE-LABEL-INDEX
          EX: twist_B -> This would be the the design group 'xvars' from the file indexed B
    """
    if historyInfo:
        # If historyInfo isn't null, then it is ready to be parsed to update the function group dropdown options
        historyInfo = json.loads(historyInfo)
        options = [{"label": i, "value": i} for i in historyInfo["optNames"]]
        return options
    else:
        # If historyInfo is null, keep the options as empty
        return []


@app.callback(
    dash.dependencies.Output("dvarChild", "options"),
    [
        dash.dependencies.Input("dvarGroup", "value"),
        dash.dependencies.Input("dvarGroup", "options"),
    ],
)
def update_dvar_child(dvarGroup, options):
    """
    This populates the 'Design Variable' dropdown based on the design groups from the 'Design Group' selected.

    Parameters
    ----------
    dvarGroup: list of str
       The selected DV groups from the 'Design Group' dropdown.

    options: list of str
        The 'Design Group' dropdown options. This is not explicitly used in this function, but is supplied
        as an input dependency so that this callback is triggered to re-update when 'refresh' is turned on. When
        'refresh' is on, once the dvarGroup options finally loads, the dvarChild's options will be update.

    Returns
    -------
    dictionary
       List of dictionaries containing the options for the 'Design Variable' dropdown.

    If multiple history files are being used - Formatted as: OPT/DV/FUNC-GROUPNAME_HIST-FILE-LABEL-INDEX_VAR-INDEX
            EX: xvars_A_0 -> This would be the first variable in the xvars group from file A
        One history file:
            EX: xvars_0 -> This would be the first variable in the xvars group
         *NOTE* If the xvars group only has one scalar quantity, the corresponding name will only be xvars with no _INDEX
            EX: If xvars only holds a scalar qty, its variable name will be simply 'xvars'

    """
    # Re-reads in history list data and passes to History API.
    try:
        # If AutoRefresh is on and the data has not yet loaded, this line will throw an error
        histList = [History(fileName) for fileName in histListArgs]
    except Exception:
        print("History file data is not processed yet")
        return []
    else:
        strlist = []
        if dvarGroup:
            for var in dvarGroup:
                var = str(var)
                varName = var
                hist = histList[0]
                if len(histList) > 1:
                    indexHist = var.split("_")[-1]
                    hist = histList[ord(indexHist) % 65]
                    varName = var[::-1].replace(indexHist + "_", "", 1)[::-1]
                varValues = hist.getValues(names=varName, major=False)
                num = len(varValues[varName][0])
                if num == 1:
                    strlist += [var]
                else:
                    strlist += [var + "_" + str(i) for i in range(num)]
        return [{"label": i, "value": i} for i in strlist]


@app.callback(
    dash.dependencies.Output("funcChild", "options"),
    [
        dash.dependencies.Input("funcGroup", "value"),
        dash.dependencies.Input("funcGroup", "options"),
    ],
)
def update_func_child(funcGroup, options):
    """
    This populates the 'Function Variable' dropdown based on the design groups from the 'Function Group' selected.

    Parameters
    ----------
    dvarGroup: list of str
       The selected DV groups from the 'Function Group' dropdown.

    options: list of str
        The 'Function Group' dropdown options. This is not explicitly used in this function, but is supplied
        as an input dependency so that this callback is triggered to re-update when 'refresh' is turned on. When
        'refresh' is on, once the funcGroup options finally loads, the funcChild's options will be update.

    Returns
    -------
    dictionary
        List of dictionaries containing the options for the 'Function Variable' dropdown.

    If multiple history files are being used - Formatted as: OPT/DV/FUNC-GROUPNAME_HIST-FILE-LABEL-INDEX_VAR-INDEX
            EX: xvars_A_0 -> This would be the first variable in the xvars group from file A
        One history file:
            EX: xvars_0 -> This would be the first variable in the xvars group
         *NOTE* If the xvars group only has one scalar quantity, the corresponding name will only be xvars with no _INDEX
            EX: Because 'obj' holds a scalar quantity, it will be displayed as simply 'obj' in the Variable dropdown.
    """
    # Re-reads in history list data and passes to History API
    try:
        # If AutoRefresh is on and the data has not yet loaded, this line will throw an error
        histList = [History(fileName) for fileName in histListArgs]
    except Exception:
        print("History file data is not processed yet")
        return []
    else:
        strlist = []
        if funcGroup:
            for var in funcGroup:
                var = str(var)
                hist = histList[0]
                varName = var
                if len(histList) > 1:
                    indexHist = var.split("_")[-1]
                    hist = histList[ord(indexHist) % 65]
                    varName = var[::-1].replace(indexHist + "_", "", 1)[::-1]
                varValues = hist.getValues(names=varName, major=False)
                num = len(varValues[varName][0])
                if num == 1:
                    strlist += [var]
                else:
                    strlist += [var + "_" + str(i) for i in range(num)]

        return [{"label": i, "value": i} for i in strlist]


@app.callback(
    dash.dependencies.Output("optChild", "options"),
    [
        dash.dependencies.Input("optGroup", "value"),
        dash.dependencies.Input("optGroup", "options"),
    ],
)
def update_opt_child(optGroup, options):
    """
    This populates the 'Optimization Variable' dropdown based on the design groups from the 'Optimization Group' selected.

    Parameters
    ----------
    optGroup: list of str
       The selected Optimization groups from the 'Optimization Group' dropdown.

    options: list of str
        The 'Optimization Group' dropdown options. This is not explicitly used in this function, but is supplied
        as an input dependency so that this callback is triggered to re-update when 'refresh' is turned on. When
        'refresh' is on, once the optGroup options finally loads, the optChild's options will be update.

    Returns
    -------
    dictionary
        List of dictionaries containing the options for the 'Optimization Variable' dropdown.
    If multiple history files are being used - Formatted as: OPT/DV/FUNC-GROUPNAME_HIST-FILE-LABEL-INDEX_VAR-INDEX

    EX: xvars_A_0 -> This would be the first variable in the xvars group from file A
    One history file:
        EX: xvars_0 -> This would be the first variable in the xvars group
        *NOTE* If the xvars group only has one scalar quantity, the corresponding name will only be xvars with no _INDEX
        EX: Because 'obj' holds a scalar quantity, it will be displayed as simply 'obj' in the Variable dropdown.
    """
    # Re-reads in history list data and passes to History API
    try:
        # If AutoRefresh is on and the data has not yet loaded, this line will throw an error
        histList = [History(fileName) for fileName in histListArgs]
    except Exception:
        print("History file data is not processed yet")
        return []
    else:
        strlist = []
        if optGroup:
            for var in optGroup:
                var = str(var)
                hist = histList[0]
                varName = var
                if len(histList) > 1:
                    indexHist = var.split("_")[-1]
                    hist = histList[ord(indexHist) % 65]
                    varName = var[::-1].replace(indexHist + "_", "", 1)[::-1]
                varValues = hist.getValues(names=varName, major=False)
                num = len(varValues[varName][0])
                if num == 1:
                    strlist += [var]
                else:
                    strlist += [var + "_" + str(i) for i in range(num)]

        return [{"label": i, "value": i} for i in strlist]


@app.callback(
    dash.dependencies.Output("dvarChild", "value"),
    [dash.dependencies.Input("dvarGroup", "value"), dash.dependencies.Input("dvarChild", "options")],
    [dash.dependencies.State("dvarChild", "value")],
)
def update_dvar_childAutoPopulate(dvarGroup, options, dvarChild):
    """
    This autopopulates the 'Design Variable' dropdown if the 'Design Group' selected has 10 or less values.

    Parameters
    ----------
    dvarGroup: list of str
       The selected Design groups from the 'Design Group' dropdown.

    options: list of str
        The 'Design Group' dropdown options. This is not explicitly used in this function, but is supplied
        as an input dependency so that this callback is triggered to re-update when 'refresh' is turned on. When
        'refresh' is on, once the dvarGroup options finally loads, the dvarChild's options will be update.

    dvarChild: list of str
        This is the current state of the selected options from the 'Design Variable' dropdown. This is used
        to check which design groups have already been selected and should not be re auto populated.

    Returns
    -------
    list of str
        List of the 'Design Variable' values to be auto selected.
    """
    # Re-reads in history list data and passes to History API
    try:
        # If AutoRefresh is on and the data has not yet loaded, this line will throw an error
        histList = [History(fileName) for fileName in histListArgs]
    except Exception:
        print("History file data is not processed yet")
        return []
    else:
        strlist = []
        if dvarGroup:
            for var in dvarGroup:
                var = str(var)
                vName = var
                hist = histList[0]
                if len(histList) > 1:
                    indexHist = var.split("_")[-1]
                    hist = histList[ord(indexHist) % 65]
                    vName = var[::-1].replace(indexHist + "_", "", 1)[::-1]
                varValues = hist.getValues(names=vName, major=False)
                # This checks if a specific group already has variables selected and shouldn't be autopopulated
                varAlreadyExists = False
                for varName in dvarChild:
                    if var in varName:
                        strlist += [varName]
                        varAlreadyExists = True
                if not varAlreadyExists:
                    if len(varValues[vName][0]) < 11:
                        num = len(varValues[vName][0])
                        if num == 1:
                            strlist += [var]
                        else:
                            strlist += [var + "_" + str(i) for i in range(num)]
        return [i for i in strlist]


@app.callback(
    dash.dependencies.Output("funcChild", "value"),
    [dash.dependencies.Input("funcGroup", "value"), dash.dependencies.Input("funcChild", "options")],
    [dash.dependencies.State("funcChild", "value")],
)
def update_func_childAutoPopulate(funcGroup, options, funcChild):
    """
    This autopopulates the 'Function Variable' dropdown if the 'Function Group' selected has 10 or less values.

    Parameters
    ----------
    funcGroup: list of str
       The selected Design groups from the 'Design Group' dropdown.

    options: list of str
        The 'Function Group' dropdown options. This is not explicitly used in this function, but is supplied
        as an input dependency so that this callback is triggered to re-update when 'refresh' is turned on. When
        'refresh' is on, once the Function options finally loads, the Function's options will be updated.

    funcChild: list of str
        This is the current state of the selected options from the 'Function Variable' dropdown. This is used
        to check which design groups have already been selected and should not be re auto populated.

    Returns
    -------
    list of str
        List of the 'Function Variable' values to be auto selected.
    """
    # Re-reads in history list data and passes to History API
    try:
        # If AutoRefresh is on and the data has not yet loaded, this line will throw an error
        histList = [History(fileName) for fileName in histListArgs]
    except Exception:
        print("History file data is not processed yet")
        return []
    else:
        strlist = []
        if funcGroup:
            for var in funcGroup:
                var = str(var)
                vName = var
                hist = histList[0]
                if len(histList) > 1:
                    indexHist = var.split("_")[-1]
                    hist = histList[ord(indexHist) % 65]
                    vName = var[::-1].replace(indexHist + "_", "", 1)[::-1]
                varValues = hist.getValues(names=vName, major=False)
                # This checks if a specific group already has variables selected and shouldn't be autopopulated
                varAlreadyExists = False
                for varName in funcChild:
                    if var in varName:
                        strlist += [varName]
                        varAlreadyExists = True
                if not varAlreadyExists:
                    if len(varValues[vName][0]) < 11:
                        num = len(varValues[vName][0])
                        if num == 1:
                            strlist += [var]
                        else:
                            strlist += [var + "_" + str(i) for i in range(num)]
        return [i for i in strlist]


@app.callback(
    dash.dependencies.Output("optChild", "value"),
    [dash.dependencies.Input("optGroup", "value"), dash.dependencies.Input("optChild", "options")],
    [dash.dependencies.State("optChild", "value")],
)
def update_opt_childAutoPopulate(optGroup, options, optChild):
    """
    This autopopulates the 'Function Variable' dropdown if the 'Optimization Group' selected has 10 or less values.

    Parameters
    ----------
    optGroup: list of str
       The selected Design groups from the 'Optimization Group' dropdown.

    options: list of str
        The 'Optimization Group' dropdown options. This is not explicitly used in this function, but is supplied
        as an input dependency so that this callback is triggered to re-update when 'refresh' is turned on. When
        'refresh' is on, once the Optimization options finally loads, the Optimization's options will be updated.

    optChild: list of str
        This is the current state of the selected options from the 'Optimization Variable' dropdown. This is used
        to check which optimization groups have already been selected and should not be re auto populated.

    Returns
    -------
    list of str
        List of the 'Optimization Variable' values to be auto selected.
    """
    # Re-reads in history list data and passes to History API
    try:
        # If AutoRefresh is on and the data has not yet loaded, this line will throw an error
        histList = [History(fileName) for fileName in histListArgs]
    except Exception:
        print("History file data is not processed yet")
        return []
    else:
        strlist = []
        if optGroup:
            for var in optGroup:
                var = str(var)
                vName = var
                hist = histList[0]
                if len(histList) > 1:
                    indexHist = var.split("_")[-1]
                    hist = histList[ord(indexHist) % 65]
                    vName = var[::-1].replace(indexHist + "_", "", 1)[::-1]
                varValues = hist.getValues(names=vName, major=False)
                # This checks if a specific group already has variables selected and shouldn't be autopopulated
                varAlreadyExists = False
                for varName in optChild:
                    if var in varName:
                        strlist += [varName]
                        varAlreadyExists = True
                if not varAlreadyExists:
                    if len(varValues[vName][0]) < 11:
                        num = len(varValues[vName][0])
                        if num == 1:
                            strlist += [var]
                        else:
                            strlist += [var + "_" + str(i) for i in range(num)]
        return [i for i in strlist]


@app.callback(
    dash.dependencies.Output("plot", "figure"),
    [
        dash.dependencies.Input("dvarGroup", "value"),
        dash.dependencies.Input("funcGroup", "value"),
        dash.dependencies.Input("optGroup", "value"),
        dash.dependencies.Input("dvarChild", "value"),
        dash.dependencies.Input("funcChild", "value"),
        dash.dependencies.Input("optChild", "value"),
        dash.dependencies.Input("plot-type", "value"),
        dash.dependencies.Input("data-type", "value"),
        dash.dependencies.Input("hidden-div", "children"),
    ],
)
def update_plot(dvarGroup, funcGroup, optGroup, dvarChild, funcChild, optChild, plotType, dataType, hiddenDiv):
    """
     This will update the plot figure accordingly in the layout if the
    design groups, function groups, design variables, function variables, plot type, or data type values are changed.
    If auto refresh is selected, it will also update the plot every time hidden-dv's values are changed, which occurs on interval.

    Parameters
    ----------
    dvarGroup : list of str
       The selected design variable groups from the DV Group dropdown.

    funcGroup: list of str
        The selected function groups (or objective) from the Function Group dropdown.

    optGroup: list of str
        The selected optimization groups from the Optimization Group dropdown.

    dvarChild: list of str
        The selected design variables from the Design Variable dropdown.

    funcChild: list of str
        The selected function variables (or objective) from the Function Variable dropdown.

    optChild: list of str
        The selected optimization variable from the Optimization Variable dropdown.

    plotType: str
        The selected plot type. (Options are currently 'Stacked' or 'Shared')

    dataType : list of str
        Contains dataType str values selected by user (i.e scale, major, delta, bounds)

    hiddenDiv: Object with the dcc.Interval component information
       This allows the plot to be updated once, or on interval if 'refresh' is selected.

    Returns
    -------
    plotly figure {} object
        The updated figure object based on new input values.
    """
    # HiddenDiv is updated once the interval component runs, so this checks if it has been updated yet and
    # the information is ready to be loaded.
    if hiddenDiv:
        # A list of the History objects for each history file argument
        histList = [History(fileName) for fileName in histListArgs]
        # List of traces to be plotted on fig
        trace = []
        fig = {}
        fig["layout"] = {}
        # List of all func groups + DV groups selected names -> used for stacked plots axis
        varsNames = []

        # Create subplots in 'fig' if stacked plot type is chosen with 1 subplot per group selected
        if plotType == "stacked":
            # Determine # of plots
            numRows = 0
            if not dvarGroup and not funcGroup and not optGroup:
                numRows = 1
            if funcGroup:
                numRows += len(funcGroup)
            if dvarGroup:
                numRows += len(dvarGroup)
            if optGroup:
                numRows += len(optGroup)
            # Determine all variable group names selected for stacked plot axis'
            varsNames = []
            dvarGroup = dvarGroup if dvarGroup else []
            funcGroup = funcGroup if funcGroup else []
            optGroup = optGroup if optGroup else []
            varsNames = dvarGroup + funcGroup + optGroup
            fig = subplots.make_subplots(
                rows=numRows, subplot_titles=varsNames, vertical_spacing=0.15, x_title="Iterations"
            )
            for i in fig["layout"]["annotations"]:
                i["font"] = dict(size=12)
            fig.update_yaxes(type=("log" if (dataType and ("log" in dataType)) else "linear"))

        # Add overall min + max traces for each DV group
        if dataType and "minMax" in dataType and dvarGroup:
            for var in dvarGroup:
                addMaxMinTraces(var, trace, dataType, histList)
                # Add min & max traces to current var's subplot if 'stacked' plot type
                if plotType == "stacked":
                    fig.append_trace(trace[len(trace) - 1], varsNames.index(var) + 1, 1)
                    fig.append_trace(trace[len(trace) - 2], varsNames.index(var) + 1, 1)

        # Add traces for each DV selected, including lower + upper bound traces if 'bounds' is selected
        if dvarChild:
            for var in dvarChild:
                groupType = "dvar"
                indexVar = var.split("_")[-1]
                varGroupName = var if var in dvarGroup else var[::-1].replace(indexVar + "_", "", 1)[::-1]
                addVarTraces(var, trace, dataType, dvarGroup, groupType)
                # Add traces to current var's subplot if 'stacked' plot type
                if plotType == "stacked":
                    if dataType and "bounds" in dataType:
                        # If bounds was selected, addVarTraces added three traces (var's trace, var's upper bound trace, var's lower bound trace)
                        fig.append_trace(trace[len(trace) - 3], varsNames.index(varGroupName) + 1, 1)
                        fig.append_trace(trace[len(trace) - 2], varsNames.index(varGroupName) + 1, 1)
                        fig.append_trace(trace[len(trace) - 1], varsNames.index(varGroupName) + 1, 1)
                    else:
                        fig.append_trace(trace[len(trace) - 1], varsNames.index(varGroupName) + 1, 1)

        # Add overall min + max traces for each function group
        if dataType and "minMax" in dataType and funcGroup:
            for var in funcGroup:
                addMaxMinTraces(var, trace, dataType, histList)
                # Add min & max traces to current var's subplot if 'stacked' plot type
                if plotType == "stacked":
                    fig.append_trace(trace[len(trace) - 1], varsNames.index(var) + 1, 1)
                    fig.append_trace(trace[len(trace) - 2], varsNames.index(var) + 1, 1)

        # Add traces for each func selected, including lower + upper bound traces if 'bounds' is selected
        if funcChild:
            for var in funcChild:
                groupType = "func"
                indexVar = var.split("_")[-1]
                varGroupName = var if var in funcGroup else var[::-1].replace(indexVar + "_", "", 1)[::-1]
                addVarTraces(var, trace, dataType, funcGroup, groupType)
                if plotType == "stacked":
                    if dataType and "bounds" in dataType:
                        fig.append_trace(trace[len(trace) - 3], varsNames.index(varGroupName) + 1, 1)
                        fig.append_trace(trace[len(trace) - 2], varsNames.index(varGroupName) + 1, 1)
                        fig.append_trace(trace[len(trace) - 1], varsNames.index(varGroupName) + 1, 1)
                    else:
                        fig.append_trace(trace[len(trace) - 1], varsNames.index(varGroupName) + 1, 1)

        # Add overall min + max traces for each optimization group
        if dataType and "minMax" in dataType and optGroup:
            for var in optGroup:
                addMaxMinTraces(var, trace, dataType, histList)
                # Add min & max traces to current var's subplot if 'stacked' plot type
                if plotType == "stacked":
                    fig.append_trace(trace[len(trace) - 1], varsNames.index(var) + 1, 1)
                    fig.append_trace(trace[len(trace) - 2], varsNames.index(var) + 1, 1)

        # Add traces for each func selected, no bounds for opt group
        if optChild:
            for var in optChild:
                groupType = "opt"
                indexVar = var.split("_")[-1]
                varGroupName = var if var in optGroup else var[::-1].replace(indexVar + "_", "", 1)[::-1]
                addVarTraces(var, trace, dataType, optGroup, groupType)
                if plotType == "stacked":
                    fig.append_trace(trace[len(trace) - 1], varsNames.index(varGroupName) + 1, 1)

        # For 'shared' plotType, set fig's 'data' key as the trace[] list, containing all the traces based on the input
        if plotType == "shared":
            fig["data"] = trace

        # Layout styling
        fig["layout"].update(
            # autosize = True,
            xaxis={
                "title": {
                    "text": "Iterations" if (plotType == "shared") else None,
                    "font": {"family": "Arial, Helvetica, sans-serif"},
                },
            },
            yaxis={
                "title": {
                    "text": None,
                },
                "type": "log" if (dataType and ("log" in dataType)) else "linear",
            },
            showlegend=True,
            font={"size": 12},
        )

        return fig

    else:
        return {}


def main():
    app.run_server(debug=True)


# Run if file is used directly, and not imported
if __name__ == "__main__":
    main()
