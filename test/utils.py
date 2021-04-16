# Standard Python modules
import unittest

# External modules
import numpy as np
from numpy.testing import assert_allclose

# First party modules
from pyoptsparse import OPT
from pyoptsparse.pyOpt_error import Error

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


# these are the informs if optimization completed successfully
SUCCESS_INFORM = {
    "SNOPT": 1,
    "IPOPT": 0,
    "NLPQLP": 0,
    "SLSQP": 0,
    "PSQP": 4,
}
# these are the name(s) of options that control the output file name
OUTPUT_FILENAMES = {
    "SNOPT": {"Print file": "_print.out", "Summary file": "_summary.out"},
    "IPOPT": {"output_file": ".out"},
    "SLSQP": {"IFILE": ".out"},
    "PSQP": {"IFILE": ".out"},
    "CONMIN": {"IFILE": ".out"},
    "NLPQLP": {"iFile": ".out"},
    "ParOpt": {"output_file": ".out"},
}


class OptTest(unittest.TestCase):
    def assert_solution(self, sol, tol):
        assert_allclose(sol.fStar, self.fStar, atol=tol, rtol=tol)
        assert_allclose(sol.xStar["xvars"], self.xStar, atol=tol, rtol=tol)
        if hasattr(sol, "lambdaStar"):
            assert_allclose(sol.lambdaStar["con"], self.lambdaStar, atol=tol, rtol=tol)

    def assert_inform(self, sol, optInform=None):
        if optInform is not None:
            self.assertEqual(sol.optInform["value"], optInform)
        else:
            # some optimizers do not have informs
            if self.optName in SUCCESS_INFORM:
                self.assertEqual(sol.optInform["value"], SUCCESS_INFORM[self.optName])

    def updateOptOptions(self, optOptions):
        optionName = OUTPUT_FILENAMES[self.optName]
        for optionName, suffix in OUTPUT_FILENAMES[self.optName].items():
            optOptions[optionName] = self._testMethodName + suffix
        return optOptions

    def optimize(self, optProb, setDV=None, optOptions=None, storeHistory=False):
        if optOptions is None:
            optOptions = {}
        # always update the output file name
        optOptions = self.updateOptOptions(optOptions)
        # Optimizer
        try:
            opt = OPT(self.optName, options=optOptions)
        except Error:
            raise unittest.SkipTest("Optimizer not available:", self.optName)

        if isinstance(setDV, str):
            optProb.setDVsFromHistory(setDV)
        elif isinstance(setDV, dict):
            optProb.setDVs(setDV)
            outDV = optProb.getDVs()
            assert_allclose(setDV["xvars"], outDV["xvars"])

        sol = opt(optProb, storeHistory=storeHistory)
        return sol
