import importlib
import inspect
import unittest
import os

class TestNoMPI(unittest.TestCase):

    # Next 3 tests check how the environment variable affects importing MPI
    def test_require_mpi(self):
        os.environ['PYOPTSPARSE_REQUIRE_MPI'] = '1'
        import pyoptsparse.pyOpt_MPI
        importlib.reload(pyoptsparse.pyOpt_MPI)
        self.assertTrue(inspect.ismodule(pyoptsparse.pyOpt_MPI.MPI))

    def test_no_mpi_requirement_given(self):
        os.environ.pop('PYOPTSPARSE_REQUIRE_MPI', None)
        import pyoptsparse.pyOpt_MPI
        importlib.reload(pyoptsparse.pyOpt_MPI)
        self.assertTrue(inspect.ismodule(pyoptsparse.pyOpt_MPI.MPI))

    def test_do_not_use_mpi(self):
        os.environ['PYOPTSPARSE_REQUIRE_MPI'] = '0'
        import pyoptsparse.pyOpt_MPI
        importlib.reload(pyoptsparse.pyOpt_MPI)
        self.assertFalse(inspect.ismodule(pyoptsparse.pyOpt_MPI.MPI))

    # Next 3 tests check how the environment variable affects using ParOpt
    def test_require_mpi_check_paropt(self):
        os.environ['PYOPTSPARSE_REQUIRE_MPI'] = '1'
        import pyoptsparse.pyParOpt.ParOpt
        importlib.reload(pyoptsparse.pyParOpt.ParOpt)
        self.assertIsNotNone(pyoptsparse.pyParOpt.ParOpt._ParOpt)

    def test_no_mpi_requirement_given_check_paropt(self):
        os.environ.pop('PYOPTSPARSE_REQUIRE_MPI', None)
        import pyoptsparse.pyParOpt.ParOpt
        importlib.reload(pyoptsparse.pyParOpt.ParOpt)
        self.assertIsNotNone(pyoptsparse.pyParOpt.ParOpt._ParOpt)

    def test_do_not_use_mpi_check_paropt(self):
        os.environ['PYOPTSPARSE_REQUIRE_MPI'] = '0'
        import pyoptsparse.pyParOpt.ParOpt
        importlib.reload(pyoptsparse.pyParOpt.ParOpt)
        self.assertIsNone(pyoptsparse.pyParOpt.ParOpt._ParOpt)

if __name__ == "__main__":
    unittest.main()


