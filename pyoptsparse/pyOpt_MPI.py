"""
pyOptSparse_MPI

A simple wrapper to MPI that enables pyOptSparse to work without
mpi4py. Only the method from the COMM object that are actually used in
pyOptSparse are included here.
"""

# Standard Python modules
import os
import warnings

# isort: off


class COMM(object):
    def __init__(self):
        self.rank = 0
        self.size = 1

    def bcast(self, obj=None, root=0):
        return obj

    def Reduce(self, sendbuf, recvbuf, op, root=0):
        recvbuf = sendbuf.copy()
        return recvbuf

    def allreduce(self, sendobj=None, recvobj=None, op=None):
        return sendobj

    def gather(self, sendobj, recvobj=None, root=0):
        return [sendobj]

    def recv(self, obj=None, source=0, tag=0, status=None):
        return obj

    def Barrier(self):
        return


class myMPI(object):
    def __init__(self):
        self.COMM_WORLD = COMM()
        self.SUM = "SUM"
        self.LOR = "OR"


# Attempt to import mpi4py.
# If PYOPTSPARSE_REQUIRE_MPI is set to a recognized positive value, attempt import
# and raise exception on failure. If set to anything else, no import is attempted.
if "PYOPTSPARSE_REQUIRE_MPI" in os.environ:
    if os.environ["PYOPTSPARSE_REQUIRE_MPI"].lower() in ["always", "1", "true", "yes"]:
        from mpi4py import MPI
    else:
        MPI = myMPI()
# If PYOPTSPARSE_REQUIRE_MPI is unset, attempt to import mpi4py, but continue on failure
# with a notification.
else:
    try:
        from mpi4py import MPI
    except ImportError:
        warn = (
            "mpi4py could not be imported. mpi4py is required to use "
            + "the parallel gradient analysis and parallel objective analysis for "
            + "non-gradient based optimizers. Continuing using a dummy MPI module "
            + "from pyOptSparse."
        )
        warnings.warn(warn)
        MPI = myMPI()
