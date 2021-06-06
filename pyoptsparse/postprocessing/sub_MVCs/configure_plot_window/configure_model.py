# --- Python 3.8 ---
"""
Model for the configure plot window
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


class ConfigurePlotModel(object):
    def __init__(self):
        self.x_var_mapping = []
        self.y_var_mapping = []

    def add_var(self, file_idx: int, axis: str):
        """
        Determine the relative index of the variable based on the file index

        Parameters
        ----------
        file_idx : int
            File index
        axis : str
            Plot axis (x or y)
        """
        count = 0
        if axis == "x":
            for i in self.x_var_mapping:
                if i[0] == file_idx:
                    count += 1

            var_idx_rel = count

            self.x_var_mapping.append([file_idx, var_idx_rel])

        elif axis == "y":
            for i in self.y_var_mapping:
                if i[0] == file_idx:
                    count += 1

            var_idx_rel = count

            self.y_var_mapping.append([file_idx, var_idx_rel])

    def remove_var(self, var_idx_abs: int, axis: str):
        if axis == "x":
            var = self.x_var_mapping.pop(var_idx_abs)

            for i in self.x_var_mapping:
                if i[0] == var[0] and i[1] > var[1]:
                    i[1] -= 1

        elif axis == "y":
            var = self.y_var_mapping.pop(var_idx_abs)

            for i in self.y_var_mapping:
                if i[0] == var[0] and i[1] > var[1]:
                    i[1] -= 1

        return var[0], var[1]
