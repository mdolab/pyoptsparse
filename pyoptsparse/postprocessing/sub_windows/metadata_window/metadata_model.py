# First party modules
from pyoptsparse.postprocessing.baseclasses.model import Model


class MetadataModel(Model):
    def __init__(self):
        """
        The metadata window model.
        """
        super(MetadataModel, self).__init__()
