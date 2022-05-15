# First party modules
from pyoptsparse.postprocessing.utils.base_classes import Model


class ConfigureModel(Model):
    def __init__(self):
        """
        Model for the configure view and controller.
        """
        super(ConfigureModel, self).__init__()
