# Standard Python modules
import importlib
import inspect
import os
import unittest

# isort: off

try:
    HAS_MPI = True
    import mpi4py  # noqa:F401
except ImportError:
    HAS_MPI = False


class TestRequireMPIEnvVar(unittest.TestCase):
    # Check how the environment variable affects importing MPI
    def test_require_mpi(self):
        if not HAS_MPI:
            raise unittest.SkipTest("mpi4py not available, skipping test.")
        os.environ["PYOPTSPARSE_REQUIRE_MPI"] = "1"
        import pyoptsparse.pyOpt_MPI

        importlib.reload(pyoptsparse.pyOpt_MPI)
        self.assertTrue(inspect.ismodule(pyoptsparse.pyOpt_MPI.MPI))

    def test_no_mpi_requirement_given(self):
        os.environ.pop("PYOPTSPARSE_REQUIRE_MPI", None)
        import pyoptsparse.pyOpt_MPI

        importlib.reload(pyoptsparse.pyOpt_MPI)
        if HAS_MPI:
            self.assertTrue(inspect.ismodule(pyoptsparse.pyOpt_MPI.MPI))
        else:
            self.assertFalse(inspect.ismodule(pyoptsparse.pyOpt_MPI.MPI))

    def test_do_not_use_mpi(self):
        os.environ["PYOPTSPARSE_REQUIRE_MPI"] = "0"
        import pyoptsparse.pyOpt_MPI

        importlib.reload(pyoptsparse.pyOpt_MPI)
        self.assertFalse(inspect.ismodule(pyoptsparse.pyOpt_MPI.MPI))


class TestRequireMPIEnvVarOnParOpt(unittest.TestCase):
    # Check how the environment variable affects using ParOpt
    def setUp(self):
        # Skip these tests if ParOpt is not installed
        try:
            from paropt import ParOpt
        except ImportError:
            raise unittest.SkipTest("Optimizer not available: paropt")

    def test_require_mpi_check_paropt(self):
        """If we require MPI, we should be able to import and instantiate a ParOpt instance."""
        os.environ["PYOPTSPARSE_REQUIRE_MPI"] = "1"
        import pyoptsparse.pyParOpt.ParOpt

        importlib.reload(pyoptsparse.pyParOpt.ParOpt)
        pyoptsparse.pyParOpt.ParOpt.ParOpt()

    def test_no_mpi_requirement_given_check_paropt(self):
        """If the environment variable is not set, we should be able to import ParOpt if we have mpi4py."""
        os.environ.pop("PYOPTSPARSE_REQUIRE_MPI", None)
        import pyoptsparse.pyParOpt.ParOpt

        importlib.reload(pyoptsparse.pyParOpt.ParOpt)
        pyoptsparse.pyParOpt.ParOpt.ParOpt()

    def test_do_not_use_mpi_check_paropt(self):
        """If we explicitly require pyoptsparse not to use MPI then we should get an import error when when try to instantiate ParOpt."""
        os.environ["PYOPTSPARSE_REQUIRE_MPI"] = "0"
        import pyoptsparse.pyParOpt.ParOpt

        importlib.reload(pyoptsparse.pyParOpt.ParOpt)
        self.assertRaises(ImportError, pyoptsparse.pyParOpt.ParOpt.ParOpt)


if __name__ == "__main__":
    unittest.main()
