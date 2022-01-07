# Standard Python modules
import unittest

# External modules
import numpy as np
from numpy.testing import assert_allclose
from pkg_resources import parse_version

# First party modules
from pyoptsparse import OPT, History
from pyoptsparse.pyOpt_error import Error

DEFAULT_TOL = 1e-12


def assertEqual(a, b):
    """
    This replaces the assertEqual functionality provided by unittest.TestCase
    so that we can call this outside the class
    """
    if not a == b:
        raise AssertionError(f"{a} and {b} are not equal.")


def assert_dict_allclose(actual, desired, atol=DEFAULT_TOL, rtol=DEFAULT_TOL, partial=False):
    """
    Simple assert for two flat dictionaries, where the values are
    assumed to be numpy arrays

    The keys are checked first to make sure that they match

    Parameters
    ----------
    actual : 1D dictionary
        A dictionary of str: array_like pairs
    desired : 1D dictionary
        A dictionary of str: array_like pairs
    atol : float
        Absolute tolerance
    rtol : float
        Relative tolerance
    partial : bool
        Whether the assertion allows for partial keys, i.e. when the dictionaries
        do not have matching keys.
        If True, we assume that ``desired`` contains a subset of the keys of ``actual``.

    Raises
    ------
    AssertionError
        if the values for any key do not match to the target tolerance.
    """
    if not partial:
        assert set(actual.keys()) == set(desired.keys())
    else:
        assert set(desired.keys()).issubset(set(actual.keys()))
    commonKeys = sorted(set(actual.keys()).intersection(set(desired.keys())))

    for key in commonKeys:
        assert_allclose(actual[key], desired[key], atol=atol, rtol=rtol, err_msg=f"Failed for key {key}")


def assert_dict_not_allclose(actual, desired, atol=DEFAULT_TOL, rtol=DEFAULT_TOL):
    """
    The opposite of assert_dict_allclose
    """
    assertEqual(set(actual.keys()), set(desired.keys()))
    for key in actual.keys():
        if np.allclose(actual[key], desired[key], atol=atol, rtol=rtol):
            raise AssertionError(f"Dictionaries are close! Got {actual} and {desired} for key {key}")


def assert_not_allclose(actual, desired, atol=DEFAULT_TOL, rtol=DEFAULT_TOL):
    """
    The numpy array version
    """
    if np.allclose(actual, desired, atol=atol, rtol=rtol):
        raise AssertionError(f"Arrays are close! Inputs are {actual} and {desired}")


def assert_optProb_size(optProb, nObj, nDV, nCon):
    """Checks that nObj, nDV and nCon are correct for optProb"""
    optProb.finalize()
    assertEqual(optProb.nObj, nObj)
    assertEqual(optProb.nCon, nCon)
    assertEqual(optProb.ndvs, nDV)


def get_dict_distance(d, d2):
    """This converts a flat 1D dict to array using sorted keys, then computes the L2 distance"""
    assertEqual(set(d.keys()), set(d2.keys()))
    a = []
    a2 = []
    for k in sorted(d.keys()):
        a.append(d[k])
        a2.append(d2[k])

    a = np.array(a)
    if a.ndim > 1:
        a = np.concatenate(a)
    a2 = np.array(a2)
    if a2.ndim > 1:
        a2 = np.concatenate(a2)
    return np.linalg.norm(a - a2)


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
    "ALPSO": {"filename": ".out"},
    "NSGA2": {},
}

# these are optimizers which are installed by default
DEFAULT_OPTIMIZERS = {"SLSQP", "PSQP", "CONMIN", "ALPSO", "NSGA2"}


class OptTest(unittest.TestCase):
    # these must be defined prior to using this class
    # name
    # optName
    # optProb
    # objfunc
    # objs
    # cons
    # DVs
    # xStar
    # fStar

    # these are optional attributes
    # sens
    # lambdaStar
    # extras

    def assert_solution_allclose(self, sol, tol, partial_x=False):
        """
        An assertion method to check that the solution object matches the expected
        optimum values defined in the class.

        Parameters
        ----------
        sol : Solution object
            The solution object after running an optimization
        tol : float
            The tolerance for assertions
        partial_x : bool, optional
            Whether partial assertion of the design variables ``x`` is allowed.
            For large problems, we may not have the full x vector available at the optimum,
            so we only check a few entries.
        """
        if not isinstance(self.xStar, list):
            self.xStar = [self.xStar]
        if not isinstance(self.fStar, list):
            self.fStar = [self.fStar]
        has_lambdaStar = hasattr(self, "lambdaStar") and hasattr(sol, "lambdaStar")
        if has_lambdaStar and not isinstance(self.lambdaStar, list):
            self.lambdaStar = [self.lambdaStar]

        if not partial_x:
            dist = []
            for x in self.xStar:
                dist.append(get_dict_distance(x, sol.xStar))
            self.sol_index = dist.index(min(dist))
        else:
            # assume we have a single solution
            self.sol_index = 0
        # now we assert against the closest solution
        # objective
        # sol.fStar was broken for earlier versions of SNOPT
        if self.optName == "SNOPT" and parse_version(self.optVersion) >= parse_version("7.7.7"):
            assert_allclose(sol.fStar, self.fStar[self.sol_index], atol=tol, rtol=tol)
        sol_objectives = np.array([sol.objectives[key].value for key in sol.objectives])
        assert_allclose(sol_objectives, self.fStar[self.sol_index], atol=tol, rtol=tol)
        # x
        assert_dict_allclose(sol.xStar, self.xStar[self.sol_index], atol=tol, rtol=tol, partial=partial_x)
        dv = sol.getDVs()
        assert_dict_allclose(dv, self.xStar[self.sol_index], atol=tol, rtol=tol, partial=partial_x)
        # lambda
        if (
            hasattr(self, "lambdaStar")
            and hasattr(sol, "lambdaStar")
            and self.lambdaStar is not None
            and sol.lambdaStar is not None
        ):
            assert_dict_allclose(sol.lambdaStar, self.lambdaStar[self.sol_index], atol=tol, rtol=tol)

    def assert_inform_equal(self, sol, optInform=None):
        """
        Check that the optInform stored in the Solution object is as expected.

        Parameters
        ----------
        sol : Solution object
            The solution object after optimization
        optInform : int, optional
            The expected inform. If None, the default inform is used, which corresponds to
            a successful optimization termination.
        """
        if optInform is not None:
            self.assertEqual(sol.optInform["value"], optInform)
        else:
            # some optimizers do not have informs
            if self.optName in SUCCESS_INFORM:
                self.assertEqual(sol.optInform["value"], SUCCESS_INFORM[self.optName])

    def get_default_hst_name(self):
        # self.id() is provided by unittest.TestCase automatically
        # which is an unambiguous identifier
        return f"{self.id()}.hst"

    def update_OptOptions_output(self, optOptions):
        """
        A general helper function to set a default output file name in optOptions
        based on the optimizer and test method name.

        Parameters
        ----------
        optOptions : dict
            The optimizer options dictionary

        Returns
        -------
        optOptions: dict
            The same optOptions dictionary, but with the output file names set.
        """
        optionName = OUTPUT_FILENAMES[self.optName]
        for optionName, suffix in OUTPUT_FILENAMES[self.optName].items():
            # self.id() is provided by unittest.TestCase automatically
            # which is an unambiguous identifier
            optOptions[optionName] = f"{self.id()}{suffix}"
        return optOptions

    def optimize(self, sens=None, setDV=None, optOptions=None, storeHistory=False, hotStart=None):
        """
        This is the workhorse routine, which will optimize self.optProb using self.optName as the optimizer.

        Parameters
        ----------
        sens : str or callable, optional
            The sens parameter which is passed to the optimizer
        setDV : str or dict, optional
            Initial DVs to use. If str, a history file path is assumed and the DVs are read from the last iteration.
            If dict, it is assumed to contain the initial DVs. None means no DVs are set and the initial values from the optProb
            are used.
        optOptions : dict, optional
            Optimizer options to use.
        storeHistory : bool or str, optional
            Controls whether a history file is written. If False, no history file is written.
            If True, a history file is written with a default name which is set to self.histFileName.
            If str, it is assumed to be the name of the history file, and self.histFileName is also assigned.
        hotStart : bool or str, optional
            Whether to use a hot start file. If False, no hot start is applied.
            If True, the default name is used as the hot start file (not storeHistory). If str, it is assumed to be the hot start file name.

        Returns
        -------
        sol : Solution object
            The solution object after optimization.

        """
        # reset counter
        self.nf = 0
        self.ng = 0
        if sens is None:
            if hasattr(self, "sens"):
                sens = self.sens

        if optOptions is None:
            optOptions = {}
        # always update the output file name
        optOptions = self.update_OptOptions_output(optOptions)
        # Optimizer
        try:
            opt = OPT(self.optName, options=optOptions)
            self.optVersion = opt.version
        except Error as e:
            if self.optName in DEFAULT_OPTIMIZERS:
                raise e
            else:
                raise unittest.SkipTest("Optimizer not available: ", self.optName)

        if isinstance(setDV, str):
            self.optProb.setDVsFromHistory(setDV)
        elif isinstance(setDV, dict):
            self.optProb.setDVs(setDV)
        DEFAULT_HST = self.get_default_hst_name()
        if storeHistory is True:
            storeHistory = DEFAULT_HST
        elif storeHistory is False:
            storeHistory = None
        if storeHistory is not None:
            self.histFileName = storeHistory

        if hotStart is True:
            hotStart = DEFAULT_HST
        elif hotStart is False:
            hotStart = None

        sol = opt(self.optProb, sens=sens, storeHistory=storeHistory, hotStart=hotStart)
        return sol

    def check_hist_file(self, tol):
        """
        We check the history file here along with the API
        """
        hist = History(self.histFileName, flag="r")
        # Metadata checks
        metadata = hist.getMetadata()
        self.assertEqual(metadata["optimizer"], self.optName)
        metadata_def_keys = [
            "optName",
            "optOptions",
            "nprocs",
            "startTime",
            "endTime",
            "optTime",
            "version",
            "optVersion",
        ]
        for key in metadata_def_keys:
            self.assertIn(key, metadata)
            # we test that SNOPT version is stored correctly
            if self.optName == "SNOPT" and key == "optVersion":
                self.assertNotEqual(metadata[key], None)

        hist.getOptProb()

        # Info checks
        self.assertEqual(set(hist.getDVNames()), self.DVs)
        self.assertEqual(set(hist.getConNames()), self.cons)
        self.assertEqual(set(hist.getObjNames()), self.objs)
        dvInfo = hist.getDVInfo()
        self.assertEqual(len(dvInfo), len(self.DVs))
        for var in self.DVs:
            self.assertEqual(dvInfo[var], hist.getDVInfo(key=var))
        conInfo = hist.getConInfo()
        self.assertEqual(len(conInfo), len(self.cons))
        objInfo = hist.getObjInfo()
        self.assertEqual(len(objInfo), len(self.objs))
        for obj in self.objs:
            self.assertEqual(objInfo[obj], hist.getObjInfo(key=obj))
        for key in ["lower", "upper", "scale"]:
            for dvName in self.DVs:
                self.assertIn(key, dvInfo[dvName])
            for con in self.cons:
                self.assertIn(key, conInfo[con])
        for obj in self.objs:
            self.assertIn("scale", objInfo[obj])

        # callCounter checks
        callCounters = hist.getCallCounters()
        last = hist.read("last")  # 'last' key should be present
        self.assertIn(last, callCounters)

        # iterKey checks
        iterKeys = hist.getIterKeys()
        for key in ["xuser", "fail", "isMajor"]:
            self.assertIn(key, iterKeys)

        # extraFuncsNames checks
        extraFuncsNames = hist.getExtraFuncsNames()
        if hasattr(self, "extras"):
            for key in self.extras:
                self.assertIn(key, extraFuncsNames)

        # getValues checks
        val = hist.getValues()

        # this check is only used for optimizers that guarantee '0' and 'last' contain funcs
        if self.optName in ["SNOPT", "PSQP"]:
            val = hist.getValues(callCounters=["0", "last"], stack=True)
            self.assertEqual(val["isMajor"].size, 2)
            self.assertTrue(val["isMajor"][0])  # the first callCounter must be a major iteration
            self.assertTrue(val["isMajor"][-1])  # the last callCounter must be a major iteration
            # check optimum stored in history file against xstar
            val = hist.getValues(callCounters="last", stack=False)
            for varName in self.DVs:
                assert_allclose(val[varName].flatten(), self.xStar[self.sol_index][varName], atol=tol, rtol=tol)

    def optimize_with_hotstart(self, tol, optOptions=None, x0=None):
        """
        This code will perform 4 optimizations, one real opt and three restarts.
        In this process, it will check various combinations of storeHistory and hotStart filenames.
        It will also call `check_hist_file` after the first optimization.
        """
        # we use a non-default starting point to test that the hotstart works
        # even if it does not match optProb initial values
        sol = self.optimize(storeHistory=True, optOptions=optOptions, setDV=x0)
        self.assert_solution_allclose(sol, tol)
        self.assertGreater(self.nf, 0)
        self.assertGreater(self.ng, 0)
        self.check_hist_file(tol)

        # re-optimize with hotstart
        sol = self.optimize(storeHistory=False, hotStart=True, optOptions=optOptions)
        self.assert_solution_allclose(sol, tol)
        # we should have zero actual function/gradient evaluations
        self.assertEqual(self.nf, 0)
        self.assertEqual(self.ng, 0)
        # another test with hotstart, this time with storeHistory = hotStart
        sol = self.optimize(storeHistory=True, hotStart=True, optOptions=optOptions)
        self.assert_solution_allclose(sol, tol)
        # we should have zero actual function/gradient evaluations
        self.assertEqual(self.nf, 0)
        self.assertEqual(self.ng, 0)
        # another test with hotstart, this time with a non-existing history file
        # this will perform a cold start
        self.optimize(storeHistory=True, hotStart="notexisting.hst", optOptions=optOptions)
        self.assertGreater(self.nf, 0)
        self.assertGreater(self.ng, 0)
        self.check_hist_file(tol)
        # final test with hotstart, this time with a different storeHistory
        sol = self.optimize(
            storeHistory=f"{self.id()}_new_hotstart.hst",
            hotStart=True,
            optOptions=optOptions,
        )
        self.assert_solution_allclose(sol, tol)
        # we should have zero actual function/gradient evaluations
        self.assertEqual(self.nf, 0)
        self.assertEqual(self.ng, 0)
