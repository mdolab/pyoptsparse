import numpy as np
from numpy.testing import assert_allclose

tol = 1e-12


def assertEqual(a, b):
    if not a == b:
        raise AssertionError(f"{a} and {b} are not equal.")


def assert_dict_allclose(actual, desired, atol=tol, rtol=tol):
    """
    Simple assert for two flat dictionaries, where the values are
    assumed to be numpy arrays

    The keys are checked first to make sure that they match
    """
    assertEqual(set(actual.keys()), set(desired.keys()))
    for key in actual.keys():
        assert_allclose(actual[key], desired[key], atol=atol, rtol=rtol)


def assert_dict_not_allclose(actual, desired, atol=tol, rtol=tol):
    """
    The opposite of assert_dict_allclose
    """
    assertEqual(set(actual.keys()), set(desired.keys()))
    for key in actual.keys():
        if np.allclose(actual[key], desired[key], atol=tol, rtol=tol):
            raise AssertionError("Dictionaries are close! Inputs are {} and {}".format(actual, desired))


def assert_not_allclose(actual, desired, atol=tol, rtol=tol):
    """
    The numpy array version
    """
    if np.allclose(actual, desired, atol=atol, rtol=tol):
        raise AssertionError("Arrays are close! Inputs are {} and {}".format(actual, desired))


def assert_optProb_size(optProb, nObj, nDV, nCon):
    """Checks that nObj, nDV and nCon are correct for optProb"""
    optProb.finalize()
    assertEqual(optProb.nObj, nObj)
    assertEqual(optProb.nCon, nCon)
    assertEqual(optProb.ndvs, nDV)
