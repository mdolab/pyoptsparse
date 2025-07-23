# Standard Python modules
import os
import sys
import unittest

# First party modules
from pyoptsparse import Optimizers, list_optimizers
from pyoptsparse.pyOpt_utils import try_import_compiled_module_from_path

# we have to unset this environment variable because otherwise
# the snopt module gets automatically imported, thus failing the import test below
os.environ.pop("PYOPTSPARSE_IMPORT_SNOPT_FROM", None)


class TestImportSnoptFromPath(unittest.TestCase):
    def test_nonexistent_path(self):
        # first unload `snopt` from namespace
        for key in list(sys.modules.keys()):
            if "snopt" in key:
                sys.modules.pop(key)
        with self.assertWarns(UserWarning):
            module = try_import_compiled_module_from_path("snopt", "/a/nonexistent/path", raise_warning=True)
            self.assertTrue(isinstance(module, str))

    def test_sys_path_unchanged(self):
        path = tuple(sys.path)
        try_import_compiled_module_from_path("snopt", "/some/path")
        self.assertEqual(tuple(sys.path), path)


class TestListOpt(unittest.TestCase):
    def test_list_optimizers(self):
        all_opt = list_optimizers()
        self.assertIn(Optimizers.SLSQP, all_opt)
        self.assertIn(Optimizers.CONMIN, all_opt)
        self.assertIn(Optimizers.PSQP, all_opt)
        self.assertIn(Optimizers.ALPSO, all_opt)
        self.assertIn(Optimizers.NSGA2, all_opt)
