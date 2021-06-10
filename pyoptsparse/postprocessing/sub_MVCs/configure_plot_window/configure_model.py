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
        self.var_mapping = []

    def add_var(self, file_idx: int):
        """
        Determine the relative index of the variable based on the file index

        Parameters
        ----------
        file_idx : int
            File index
        """
        count = 0
        for i in self.var_mapping:
            if i[0] == file_idx:
                count += 1

        var_idx_rel = count

        self.var_mapping.append([file_idx, var_idx_rel])

    def remove_var(self, var_idx_abs: int):

        var = self.var_mapping.pop(var_idx_abs)

        for i in self.var_mapping:
            if i[0] == var[0] and i[1] > var[1]:
                i[1] -= 1

        return var[0], var[1]
