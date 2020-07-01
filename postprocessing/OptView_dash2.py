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



# Save history file arguments to histList
major_python_version = sys.version_info[0]
parser = argparse.ArgumentParser()
parser.add_argument(
    "histFile", nargs="*", type=str, default="opt_hist.hst", help="Specify the history file to be plotted"
)
args = parser.parse_args()
histList = args.histFile



class ReadOptHist(OVBaseClass):

    def __init__(self, histList):
        # Initialize variables and save inputs from user
        self.histList = histList
#Will i be doing bolor bounds?
        self.color_bounds = [0.0, 0.0]


    # Call baseclass function that sets up and initializes needed inputs from db
    # Initailzies the following lists,dicts
    # self.func_data_all,func_major_all
    # self.var_data_all,var_major_all, and self.scaling,bounds 
    # from history file data to save major/minor iteration info
    self.OptimizationHistory()

#Look into this!
    def refresh_history(self):
        """
        Refresh opt_his data if the history file has been updated.
        """
        self.OptimizationHistory()



# ======================================================================
# Dash GUI Set-UP
# ======================================================================

# Intializes Opt object to hold optimziation history information
Opt = ReadOptHist(histList)





