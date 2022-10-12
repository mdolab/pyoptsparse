# First party modules
from pyoptsparse.postprocessing.baseclasses.model import Model


class ConfigureModel(Model):
    def __init__(self):
        """
        Model for the configure view and controller.
        """
        super(ConfigureModel, self).__init__()
