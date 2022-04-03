# Standard Python modules

# External modules

# Local modules
from pyoptsparse.postprocessing.utils.base_classes import Model


class MetadataModel(Model):
    def __init__(self):
        """
        The metadata window model.
        """
        super(MetadataModel, self).__init__()
