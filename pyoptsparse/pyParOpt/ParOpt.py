# First party modules
from pyoptsparse.pyOpt_error import Error
from pyoptsparse.pyOpt_optimizer import Optimizer

try:
    # External modules
    from paropt.paropt_pyoptsparse import ParOptSparse as ParOpt
except ImportError:

    class ParOpt(Optimizer):
        def __init__(self, raiseError=True, options={}):
            name = "ParOpt"
            category = "Local Optimizer"
            self.set_options = {}
            self.informs = {}
            super().__init__(
                name,
                category,
                defaultOptions=self.defOpts,
                informs=self.informs,
                options=options,
            )
            if raiseError:
                raise Error("There was an error importing ParOpt")
