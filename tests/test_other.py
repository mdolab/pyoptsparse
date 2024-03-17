# Standard Python modules
import sys
import unittest

# First party modules
from pyoptsparse.pyOpt_utils import try_import_compiled_module_from_path  # noqa: E402


class TestImportSnoptFromPath(unittest.TestCase):
    def test_nonexistent_path(self):
        with self.assertWarns(ImportWarning):
            self.assertIsNone(try_import_compiled_module_from_path("snopt", "/a/nonexistent/path"))

    def test_sys_path_unchanged(self):
        path = tuple(sys.path)
        try_import_compiled_module_from_path("snopt", "/some/path")
        self.assertEqual(tuple(sys.path), path)
