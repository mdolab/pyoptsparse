# Standard Python modules
import os

# First party modules
from pyoptsparse.pyOpt_optimizer import Optimizer

NAME = "ParOpt"
CATEGORY = "Local Optimizer"

# Attempt to import ParOpt/mpi4py
# If PYOPTSPARSE_REQUIRE_MPI is not set to a recognized positive value, that means the user is explicitly requiring
# pyOptSparse not to use MPI. In this case, we cannot use ParOpt because it requires mpi4py.
if "PYOPTSPARSE_REQUIRE_MPI" in os.environ and os.environ["PYOPTSPARSE_REQUIRE_MPI"].lower() not in [
    "always",
    "1",
    "true",
    "yes",
]:

    class ParOpt(Optimizer):
        def __init__(self, raiseError=True, options={}, sparse=True):
            super().__init__(
                NAME,
                CATEGORY,
            )
            if raiseError:
                raise ImportError(
                    "ParOpt was not imported, as requested by the environment variable 'PYOPTSPARSE_REQUIRE_MPI'"
                )


# In all other cases, we attempt to import ParOpt and mpi4py. If either import fails, we disable the optimizer.
else:
    try:
        # External modules
        from mpi4py import MPI  # noqa:F401

        try:
            # External modules
            from paropt.paropt_pyoptsparse import ParOptSparse as ParOpt
        except ImportError:

            class ParOpt(Optimizer):
                def __init__(self, raiseError=True, options={}, sparse=True):
                    super().__init__(
                        NAME,
                        CATEGORY,
                    )
                    if raiseError:
                        raise ImportError("There was an error importing ParOpt")

    except ImportError:

        class ParOpt(Optimizer):
            def __init__(self, raiseError=True, options={}, sparse=True):
                super().__init__(
                    NAME,
                    CATEGORY,
                )
                if raiseError:
                    raise ImportError("There was an error importing mpi4py, which is required for ParOpt")
