# First party modules
from pyoptsparse.pyOpt_optimizer import Optimizer

try:
    # External modules
    from paropt.paropt_pyoptsparse import ParOptSparse as ParOpt
except ImportError as e:

    def make_cls(e):
        class ParOpt(Optimizer):
            def __init__(self, raiseError=True, options={}):
                name = "ParOpt"
                category = "Local Optimizer"
                self.defOpts = {}
                self.informs = {}
                super().__init__(
                    name,
                    category,
                    defaultOptions=self.defOpts,
                    informs=self.informs,
                    options=options,
                )
                if raiseError:
                    raise e

        return ParOpt

    ParOpt = make_cls(e)
