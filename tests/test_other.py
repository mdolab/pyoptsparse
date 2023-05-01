# Standard Python modules
import sys
import unittest

# First party modules
from pyoptsparse.pySNOPT.pySNOPT import _import_snopt_from_path


class TestImportSnoptFromPath(unittest.TestCase):
    def test_nonexistent_path(self):
        with self.assertWarns(ImportWarning):
            self.assertIsNone(_import_snopt_from_path("/a/nonexistent/path"))

    def test_sys_path_unchanged(self):
        path = tuple(sys.path)
        _import_snopt_from_path("/some/path")
        self.assertEqual(tuple(sys.path), path)
