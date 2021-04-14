# --- Python 3.8 ---
"""
Encapsulates all state control for PyQt5 widgets
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


class StateController:
    def __init__(self, view):
        self._view = view

    def setInitialState(self):
        self._view.x_cbox.setDisabled(True)
        self._view.x_undo_btn.setDisabled(True)
        self._view.x_clear_btn.setDisabled(True)

        self._view.y_cbox.setDisabled(True)
        self._view.y_undo_btn.setDisabled(True)
        self._view.y_clear_btn.setDisabled(True)

        self._view.stack_plot_opt.setDisabled(True)
        self._view.abs_delta_opt.setDisabled(True)
        self._view.minor_itr_opt.setDisabled(True)
        self._view.major_itr_opt.setDisabled(True)
        self._view.min_max_opt.setDisabled(True)
        self._view.bound_opt.setDisabled(True)

        self._view.plot_btn.setDisabled(True)
        self._view.clear_plot_btn.setDisabled(True)
        self._view.refresh_btn.setDisabled(True)

        self._view.scale_var_togg.setDisabled(True)
        self._view.auto_refresh_togg.setDisabled(True)

    def setAddFileState(self):
        # --- State control after file loads ---
        self._view.x_cbox.setDisabled(False)
        self._view.x_undo_btn.setDisabled(False)
        self._view.x_clear_btn.setDisabled(False)

        self._view.y_cbox.setDisabled(False)
        self._view.y_undo_btn.setDisabled(False)
        self._view.y_clear_btn.setDisabled(False)

        self._view.minor_itr_opt.setDisabled(False)
        self._view.major_itr_opt.setDisabled(False)

        self._view.refresh_btn.setDisabled(False)

        self._view.auto_refresh_togg.setDisabled(False)

    def setMajorMinorIterCheckedState(self):
        self._view.x_cbox.setDisabled(True)
        self._view.x_label.setStyleSheet("background-color: rgba(192, 192, 192, 10); border: 1px solid black;")
        self._view.x_undo_btn.setDisabled(True)
        self._view.x_clear_btn.setDisabled(True)

    def setMajorMinorIterUncheckedState(self):
        self._view.x_cbox.setDisabled(False)
        self._view.x_label.setStyleSheet("background-color: white; border: 1px solid black;")
        self._view.x_undo_btn.setDisabled(False)
        self._view.x_clear_btn.setDisabled(False)

    def setStackedPlotState(self, state: bool):
        self._view.stack_plot_opt.setDisabled(state)

    def setClearVarState(self):
        self._view.stack_plot_opt.setDisabled(True)
        self._view.abs_delta_opt.setDisabled(True)
        self._view.min_max_opt.setDisabled(True)
        self._view.bound_opt.setDisabled(True)

    def setAbsDeltaState(self, state: bool):
        self._view.abs_delta_opt.setDisabled(state)

    def setPlotState(self, state: bool):
        self._view.plot_btn.setDisabled(state)
        self._view.clear_plot_btn.setDisabled(state)
        self._view.bound_opt.setDisabled(state)
        self._view.scale_var_togg.setDisabled(state)
