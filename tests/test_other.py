# Standard Python modules
import os
import sys
import unittest

# First party modules
from pyoptsparse import Optimizers, list_optimizers
from pyoptsparse.pyOpt_solution import SolutionInform
from pyoptsparse.pyOpt_utils import import_module

# we have to unset this environment variable because otherwise
# the snopt module gets automatically imported, thus failing the import test below
os.environ.pop("PYOPTSPARSE_IMPORT_SNOPT_FROM", None)


class TestImportSnoptFromPath(unittest.TestCase):
    def test_nonexistent_path(self):
        # first unload `snopt` from namespace
        for key in list(sys.modules.keys()):
            if "snopt" in key:
                sys.modules.pop(key)
        with self.assertRaises(ImportError):
            import_module("snopt", ["/a/nonexistent/path"], on_error="raise")

    def test_import_standard(self):
        loaded = import_module("os")
        assert loaded.__name__ == "os"

    def test_import_nonexistent(self):
        with self.assertRaises(ImportError):
            _ = import_module("not_a_module", on_error="raise")

        e = import_module("not_a_module", on_error="return")
        assert isinstance(e, Exception)
        assert "No module" in str(e)

    def test_sys_path_unchanged(self):
        path = tuple(sys.path)
        import_module("somemodule", ["/some/path"])
        self.assertEqual(tuple(sys.path), path)


class TestListOpt(unittest.TestCase):
    def test_list_optimizers(self):
        all_opt = list_optimizers()
        self.assertIn(Optimizers.SLSQP, all_opt)
        self.assertIn(Optimizers.CONMIN, all_opt)
        self.assertIn(Optimizers.PSQP, all_opt)
        self.assertIn(Optimizers.ALPSO, all_opt)
        self.assertIn(Optimizers.NSGA2, all_opt)


class TestSolInform(unittest.TestCase):
    def test_sol_inform_key_access(self):
        sol = SolutionInform(value=1, message="test message")
        assert sol.value == 1
        assert sol["value"] == 1
        assert sol.message == "test message"
        assert sol["text"] == "test message"
