# --- Python 3.8 ---
"""
Controller for the main view.  Interacts with the data models and
handles all user input and response functionality.  Controller can
only update the view based on user input.  If a view state is changed
which requires a messagebox view, that view is created by the controller
but managed seperately.
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


class MainController:
    """
    Contains functionality for user input and software
    response for the main view.
    """

    def __init__(self, view):
        self.view = view
