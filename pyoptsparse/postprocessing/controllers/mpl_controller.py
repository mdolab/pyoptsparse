# --- Python 3.8 ---
"""
Controller for all matplotlib canvas related input and response
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


class MplController:
    """Handles all matplotlib input and updates the view accordingly."""

    def __init__(self, canvas):
        self.canvas = canvas
