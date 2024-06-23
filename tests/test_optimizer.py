"""Test for Optimizer"""

# Standard Python modules
import unittest

# External modules
import numpy as np
from parameterized import parameterized

# First party modules
from pyoptsparse import OPT, Optimization

MASTERFUNC_OUTPUTS = ["fobj", "fcon", "gobj", "gcon"]


class TestOptimizer(unittest.TestCase):
    tol = 1e-12

    def get_objfunc(self, failFlag=False):
        """
        Return an objfunc callable function where we can choose whether
        the fail flag will be returned as True or False.
        """
        # Initialize iters to infinite so the fail flag is never thrown on setup
        iters = np.inf

        def objfunc(xdict):
            """
            This is a simple quadratic test function with linear constraints.
            The actual problem doesn't really matter, since we are not testing optimization,
            but just optProb. However, we need to initialize and run an optimization
            in order to have optimizer-specific fields in optProb populated, such as
            jacIndices.

            This problem is probably not feasible, but that's okay.

            Reset the iteration counter with a special call that includes
            a nonfinite value in the design variable vector.
            """
            funcs = {}
            funcs["obj_0"] = 0
            for x in xdict.keys():
                funcs["obj_0"] += np.sum(np.power(xdict[x], 2))
            for iCon, nc in enumerate(self.nCon):
                conName = f"con_{iCon}"
                funcs[conName] = np.zeros(nc)
                for x in xdict.keys():
                    for j in range(nc):
                        funcs[conName][j] = (iCon + 1) * np.sum(xdict[x])

            # Throw the fail flag if it's in the specified range or True
            nonlocal iters
            if isinstance(failFlag, tuple):
                if not len(failFlag) == 2:
                    raise ValueError("Fail flag must be a tuple of (iter start fail, iter end fail) or a boolean")
                fail = failFlag[0] <= iters < failFlag[1]
            elif isinstance(failFlag, bool):
                fail = failFlag
            else:
                raise ValueError("Fail flag must be a tuple of (iter start fail, iter end fail) or a boolean")
            iters += 1

            # Reset iteration counter if any non-finite values in DV dict
            for xVec in xdict.values():
                if not np.all(np.isfinite(xVec)):
                    iters = 0
                    break
            return funcs, fail

        return objfunc

    def setup_optProb(self, failFlag=False, nObj=1, nDV=[4], nCon=[2]):
        """
        This function sets up a general optimization problem, with arbitrary
        DVs, constraints and objectives.
        Arbitrary scaling for the various parameters can also be specified.
        """
        self.nObj = nObj
        self.nDV = nDV
        self.nCon = nCon

        # Optimization Object
        self.optProb = Optimization("Configurable Test Problem", self.get_objfunc(failFlag=failFlag))
        self.x0 = {}
        # Design Variables
        for iDV in range(len(nDV)):
            n = nDV[iDV]
            x0 = np.ones(n)
            dvName = f"x{iDV}"
            self.x0[dvName] = x0
            self.optProb.addVarGroup(
                dvName,
                n,
                lower=-1,
                upper=1,
                value=x0,
            )

        # Constraints
        for iCon in range(len(nCon)):
            nc = nCon[iCon]
            self.optProb.addConGroup(
                f"con_{iCon}",
                nc,
                lower=-5,
                upper=5,
            )

        # Objective
        for iObj in range(nObj):
            self.optProb.addObj(f"obj_{iObj}")

        # Finalize
        self.optProb.printSparsity()
        # create and store optimizer
        self.opt = OPT("slsqp", options={"IFILE": "optProb_SLSQP.out"})
        self.opt(self.optProb, sens="FD")

        # Call the masterFunc with some infinite DVs so it resets iters
        self.opt._masterFunc(np.full(np.sum(nDV), np.inf), ["fobj"])

    def test_masterFunc_fobj_fail(self):
        """
        Test that if the objective fails when _masterFunc is called,
        the fail flag is returned with the expected value.
        """
        nDV = [4]
        self.setup_optProb(failFlag=(1, 100), nDV=nDV)

        x = np.ones(np.sum(nDV), dtype=float)

        # Do not fail
        _, fail = self.opt._masterFunc(x, ["fobj"])
        self.assertFalse(fail)

        # Should fail on the second function call
        x += 1  # change x so it doesn't use the cache
        _, fail = self.opt._masterFunc(x, ["fobj"])
        self.assertTrue(fail)

    @parameterized.expand(MASTERFUNC_OUTPUTS)
    def test_masterFunc_output_fail_cache(self, output):
        """
        Test that if an output fails when _masterFunc is called
        and it is then called again with the same x vector,
        the fail flag is returned with the expected value.
        """
        nDV = [4]
        # Set fail flag to (0, 1) so we know for sure that it's using
        # the cache since the only failure is on the first call
        self.setup_optProb(failFlag=(0, 1), nDV=nDV)

        x = np.ones(np.sum(nDV), dtype=float)

        # Fail
        _, fail = self.opt._masterFunc(x, [output])
        self.assertTrue(fail)

        # Should fail with the same x vector using the cache
        _, fail = self.opt._masterFunc(x, [output])
        self.assertTrue(fail)

        # Do the same thing one more time to make sure the cache is really really working
        _, fail = self.opt._masterFunc(x, [output])
        self.assertTrue(fail)

    def test_masterFunc_gobj_fail_cache(self):
        """
        Test that if the gradient fails when _masterFunc is called
        and it is then called again with the same x vector,
        the fail flag is returned with the expected value.
        """
        nDV = [4]
        self.setup_optProb(failFlag=True, nDV=nDV)

        x = np.ones(np.sum(nDV), dtype=float)

        # Fail
        _, _, fail = self.opt._masterFunc(x, ["gcon", "gobj"])
        self.assertTrue(fail)

        # Should fail with the same x vector using the cache
        _, fail = self.opt._masterFunc(x, ["gobj"])
        self.assertTrue(fail)

    def test_masterFunc_fobj_fcon_cache_fail(self):
        """
        Test that if the objective fails when _masterFunc is called
        and then the constraints are called, it still returns a failure.
        """
        nDV = [4]
        self.setup_optProb(failFlag=(1, 100), nDV=nDV)

        x = np.ones(np.sum(nDV), dtype=float)

        # Do not fail
        _, fail = self.opt._masterFunc(x, ["fobj"])
        self.assertFalse(fail)

        # Check that the cached value does not fail either
        _, fail = self.opt._masterFunc(x, ["fcon"])
        self.assertFalse(fail)

        # Should fail on the second function call
        x += 1  # change x so it doesn't use the cache
        _, fail = self.opt._masterFunc(x, ["fobj"])
        self.assertTrue(fail)

        # Check that the cached value now fails too
        _, fail = self.opt._masterFunc(x, ["fcon"])
        self.assertTrue(fail)

    def test_masterFunc_fail_then_success(self):
        """
        Test that if the objective/constraint fails when _masterFunc is called
        and then it succeeds, the fail flag is no longer true.
        """
        nDV = [4, 5]
        self.setup_optProb(failFlag=(0, 1), nDV=nDV)

        x = np.ones(np.sum(nDV), dtype=float)

        # Fail
        _, _, fail = self.opt._masterFunc(x, ["fobj", "fcon"])
        self.assertTrue(fail)

        # Should succeed on the second call
        x += 1  # change x so it doesn't use the cache
        _, _, fail = self.opt._masterFunc(x, ["fobj", "fcon"])
        self.assertFalse(fail)

    def test_masterFunc_fail_grad_after_fail_func(self):
        """
        Test that if the _masterFunc is called to compute the gradients on
        an x that isn't in the cache and the primal fails, it returns a
        fail flag for the gradient too.
        """
        nDV = [4, 5]
        self.setup_optProb(failFlag=True, nDV=nDV)

        x = np.ones(np.sum(nDV), dtype=float) + 5

        # Fail obj gradient on DVs that haven't been evaluated when the primal fails
        _, fail = self.opt._masterFunc(x, ["gobj"])
        self.assertTrue(fail)

        # Fail con gradient on DVs that haven't been evaluated when the primal fails
        x += 1
        _, fail = self.opt._masterFunc(x, ["gcon"])
        self.assertTrue(fail)

    def test_masterFunc_succeed_grad_after_fail_func(self):
        """
        Test that if the _masterFunc is called to compute the gradients on
        an x that is in the cache and the primal fails, it returns a
        False fail flag for the gradient.
        """
        nDV = [4, 5]
        self.setup_optProb(failFlag=(0, 1), nDV=nDV)

        x = np.ones(np.sum(nDV), dtype=float) + 5

        _, fail = self.opt._masterFunc(x, ["fobj"])  # call primal to put x in the cache
        self.assertTrue(fail)

        # Gradient succeeds even though primal failed
        _, _, fail = self.opt._masterFunc(x, ["gobj", "gcon"])
        self.assertFalse(fail)


if __name__ == "__main__":
    unittest.main()
