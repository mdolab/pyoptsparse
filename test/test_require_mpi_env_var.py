import importlib
import inspect
import unittest
import os
import sys

if sys.version_info[0] == 2:
    reload_func = reload  # noqa: F821
else:
    reload_func = importlib.reload


class TestRequireMPIEnvVar(unittest.TestCase):
    def setUp(self):
        # check if mpi4py is installed
        try:
            import mpi4py  # noqa:F401
        except ImportError:
            raise unittest.SkipTest("No mpi4py available, skipping MPI tests.")

    # Check how the environment variable affects importing MPI
    def test_require_mpi(self):
        os.environ["PYOPTSPARSE_REQUIRE_MPI"] = "1"
        import pyoptsparse.pyOpt_MPI

        reload_func(pyoptsparse.pyOpt_MPI)
        self.assertTrue(inspect.ismodule(pyoptsparse.pyOpt_MPI.MPI))

    def test_no_mpi_requirement_given(self):
        os.environ.pop("PYOPTSPARSE_REQUIRE_MPI", None)
        import pyoptsparse.pyOpt_MPI

        reload_func(pyoptsparse.pyOpt_MPI)
        self.assertTrue(inspect.ismodule(pyoptsparse.pyOpt_MPI.MPI))

    def test_do_not_use_mpi(self):
        os.environ["PYOPTSPARSE_REQUIRE_MPI"] = "0"
        import pyoptsparse.pyOpt_MPI

        reload_func(pyoptsparse.pyOpt_MPI)
        self.assertFalse(inspect.ismodule(pyoptsparse.pyOpt_MPI.MPI))


class TestRequireMPIEnvVarOnParOpt(unittest.TestCase):
    # Check how the environment variable affects using ParOpt
    def setUp(self):
        # Just check to see if ParOpt is installed before doing any testing
        try:
            from paropt import ParOpt as _ParOpt  # noqa: F401
        except ImportError:
            raise unittest.SkipTest("Optimizer not available:", "paropt")

    def test_require_mpi_check_paropt(self):
        os.environ["PYOPTSPARSE_REQUIRE_MPI"] = "1"
        import pyoptsparse.pyParOpt.ParOpt

        reload_func(pyoptsparse.pyParOpt.ParOpt)
        self.assertIsNotNone(pyoptsparse.pyParOpt.ParOpt._ParOpt)

    def test_no_mpi_requirement_given_check_paropt(self):
        os.environ.pop("PYOPTSPARSE_REQUIRE_MPI", None)
        import pyoptsparse.pyParOpt.ParOpt

        reload_func(pyoptsparse.pyParOpt.ParOpt)
        self.assertIsNotNone(pyoptsparse.pyParOpt.ParOpt._ParOpt)

    def test_do_not_use_mpi_check_paropt(self):
        os.environ["PYOPTSPARSE_REQUIRE_MPI"] = "0"
        import pyoptsparse.pyParOpt.ParOpt

        reload_func(pyoptsparse.pyParOpt.ParOpt)
        self.assertIsNone(pyoptsparse.pyParOpt.ParOpt._ParOpt)


if __name__ == "__main__":
    unittest.main()
