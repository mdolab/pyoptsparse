# --- Python 3.8 ---
"""
Data structure and management class for history files.  The controller
only has access to top level data for plotting.  Data manipulation
only occurs here and not in the controller or views.
"""

# ==============================================================================
# Standard Python modules
# ==============================================================================

# ==============================================================================
# External Python modules
# ==============================================================================

# ==============================================================================
# Extension modules
# ==============================================================================
from pyoptsparse import History


class HistoryFileModel:
    """Manages top-level data for the controller"""

    def __init__(self, fp):
        self._file = History(fileName=fp)

    def changeFile(self, fp):
        self._file = History(fileName=fp)
