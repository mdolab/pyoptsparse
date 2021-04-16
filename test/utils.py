import unittest
import numpy as np
from numpy.testing import assert_allclose
from pyoptsparse.pyOpt_error import Error
from pyoptsparse import OPT

DEFAULT_TOL = 1e-12


def assertEqual(a, b):
    """This replaces the assertEqual functionality provided by unittest.TestCase"""
    if not a == b:
        raise AssertionError(f"{a} and {b} are not equal.")


def assert_dict_allclose(actual, desired, atol=DEFAULT_TOL, rtol=DEFAULT_TOL):
    """
    Simple assert for two flat dictionaries, where the values are
    assumed to be numpy arrays

    The keys are checked first to make sure that they match
    """
    assertEqual(set(actual.keys()), set(desired.keys()))
    for key in actual.keys():
        assert_allclose(actual[key], desired[key], atol=atol, rtol=rtol)


def assert_dict_not_allclose(actual, desired, atol=DEFAULT_TOL, rtol=DEFAULT_TOL):
    """
    The opposite of assert_dict_allclose
    """
    assertEqual(set(actual.keys()), set(desired.keys()))
    for key in actual.keys():
        if np.allclose(actual[key], desired[key], atol=atol, rtol=rtol):
            raise AssertionError("Dictionaries are close! Inputs are {} and {}".format(actual, desired))


def assert_not_allclose(actual, desired, atol=DEFAULT_TOL, rtol=DEFAULT_TOL):
    """
    The numpy array version
    """
    if np.allclose(actual, desired, atol=atol, rtol=rtol):
        raise AssertionError("Arrays are close! Inputs are {} and {}".format(actual, desired))


def assert_optProb_size(optProb, nObj, nDV, nCon):
    """Checks that nObj, nDV and nCon are correct for optProb"""
    optProb.finalize()
    assertEqual(optProb.nObj, nObj)
    assertEqual(optProb.nCon, nCon)
    assertEqual(optProb.ndvs, nDV)


class OptTest(unittest.TestCase):

    # these are the informs if optimization completed successfully
    optInform = {
        "SNOPT": 1,
        "IPOPT": 0,
        "NLPQLP": 0,
        "SLSQP": 0,
        "PSQP": 4,
    }

    def check_solution(self, sol, tol):
        assert_allclose(sol.fStar, self.fStar, atol=tol, rtol=tol)
        assert_allclose(sol.xStar["xvars"], self.xStar, atol=tol, rtol=tol)
        if hasattr(sol, "lambdaStar"):
            assert_allclose(sol.lambdaStar["con"], self.lambdaStar, atol=tol, rtol=tol)

    def check_inform(self, sol, optName, optInform=None):
        if optInform is not None:
            self.assertEqual(sol.optInform["value"], optInform)
        else:
            # some optimizers do not have informs
            if optName in self.optInform:
                self.assertEqual(sol.optInform["value"], self.optInform[optName])

    def updateOptOptions(self, optName, optOptions):
        fileNameMapping = {
            "SNOPT": {"Print file": "_print.out", "Summary file": "_summary.out"},
            "IPOPT": {"output_file": ".out"},
            "SLSQP": {"IFILE": ".out"},
            "PSQP": {"IFILE": ".out"},
            "CONMIN": {"IFILE": ".out"},
            "NLPQLP": {"iFile": ".out"},
            "ParOpt": {"output_file": ".out"},
        }
        optionName = fileNameMapping[optName]
        for optionName, suffix in fileNameMapping[optName].items():
            optOptions[optionName] = self._testMethodName + suffix
        return optOptions

    def optimize(self, optProb, optName, setDV=None, optOptions={}, storeHistory=False):
        # Optimizer
        try:
            opt = OPT(optName, options=optOptions)
        except Error:
            raise unittest.SkipTest("Optimizer not available:", optName)

        if isinstance(setDV, str):
            optProb.setDVsFromHistory(setDV)
        elif isinstance(setDV, dict):
            optProb.setDVs(setDV)
            outDV = optProb.getDVs()
            assert_allclose(setDV["xvars"], outDV["xvars"])

        sol = opt(optProb, storeHistory=storeHistory)
        return sol
