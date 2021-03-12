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


class HistFile:
    """
    Data structure and helpful functionality for loading and
    parsing pyOptSparse history files
    """

    def __init__(self, fp: str):
        """
        Initializer for the data model.

        Parameters
        ----------
        fp : str
            File path to the history file
        """
        self.fp = fp


class DataManager:
    """Manages top-level data for the controller"""

    def __init__(self):
        pass
