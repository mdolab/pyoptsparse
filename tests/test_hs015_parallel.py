"""Test solution of problem HS15 from the Hock & Schittkowski collection"""

# Standard Python modules
import unittest

# External modules
import numpy as np

try:
    HAS_MPI = True
    # External modules
    from mpi4py import MPI

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
except ImportError:
    HAS_MPI = False

# First party modules
from pyoptsparse import Optimization

# Local modules
from testing_utils import OptTest


class TestHS15(OptTest):
    ## Solve test problem HS15 from the Hock & Schittkowski collection.
    #
    #  min   100 (x2 - x1^2)^2 + (1 - x1)^2
    #  s.t.  x1 x2 >= 1
    #        x1 + x2^2 >= 0
    #        x1 <= 0.5
    #
    #  The standard start point (-2, 1) usually converges to the standard
    #  minimum at (0.5, 2.0), with final objective = 306.5.
    #  Sometimes the solver converges to another local minimum
    #  at (-0.79212, -1.26243), with final objective = 360.4.
    ##

    if not HAS_MPI:
        raise unittest.SkipTest("MPI not available")

    N_PROCS = 2  # Run case on two procs

    name = "HS015"
    DVs = {"xvars"}
    cons = {"con"}
    objs = {"obj"}
    extras = {"extra1", "extra2"}
    fStar = [
        306.5,
        360.379767,
    ]
    xStar = [
        {"xvars": (0.5, 2.0)},
        {"xvars": (-0.79212322, -1.26242985)},
    ]
    optOptions = {}

    def objfunc(self, xdict):
        self.nf += 1
        x = xdict["xvars"]
        funcs = {}
        funcs["obj"] = [100 * (x[1] - x[0] ** 2) ** 2 + (1 - x[0]) ** 2]
        conval = np.zeros(2, "D")
        conval[0] = x[0] * x[1]
        conval[1] = x[0] + x[1] ** 2
        funcs["con"] = conval
        # extra keys
        funcs["extra1"] = 0.0
        funcs["extra2"] = 1.0
        fail = False
        return funcs, fail

    def sens(self, xdict, funcs):
        self.ng += 1
        x = xdict["xvars"]
        funcsSens = {}
        funcsSens["obj"] = {
            "xvars": [2 * 100 * (x[1] - x[0] ** 2) * (-2 * x[0]) - 2 * (1 - x[0]), 2 * 100 * (x[1] - x[0] ** 2)]
        }
        funcsSens["con"] = {"xvars": [[x[1], x[0]], [1, 2 * x[1]]]}
        fail = False
        return funcsSens, fail

    def setup_optProb(self):
        # Optimization Object
        self.optProb = Optimization("HS15 Constraint Problem", self.objfunc)

        # Design Variables
        lower = [-5.0, -5.0]
        upper = [0.5, 5.0]
        value = [-2, 1.0]
        self.optProb.addVarGroup("xvars", 2, lower=lower, upper=upper, value=value)

        # Constraints
        lower = [1.0, 0.0]
        upper = [None, None]
        self.optProb.addConGroup("con", 2, lower=lower, upper=upper)

        # Objective
        self.optProb.addObj("obj")

    @staticmethod
    def my_snstop(iterDict):
        """manually terminate SNOPT after 1 major iteration if"""
        # print("Iteration", iterDict["nMajor"])
        # print("Rank", rank)
        return_idx = 0
        if rank == 1 and iterDict["nMajor"] == 1:
            return_idx = comm.bcast(1, root=1)
        return_idx = comm.bcast(return_idx, root=1)
        # comm.Barrier()
        # print(f"return_iDX on {rank}", return_idx)
        return return_idx

    def test_optimization(self):
        self.optName = "SNOPT"
        self.setup_optProb()
        sol = self.optimize()
        # Check Solution
        self.assert_solution_allclose(sol, 1e-12)
        # Check informs
        self.assert_inform_equal(sol)

    def test_snopt_snstop(self):
        self.optName = "SNOPT"
        self.setup_optProb()
        optOptions = {
            "snSTOP function handle": self.my_snstop,
        }
        sol = self.optimize(optOptions=optOptions, storeHistory=True)
        # Check informs
        # we should get 70/74
        self.assert_inform_equal(sol, optInform=74)


if __name__ == "__main__":
    unittest.main()
