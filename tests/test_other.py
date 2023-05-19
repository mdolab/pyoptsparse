# Standard Python modules
import os
import sys
import unittest

# we have to unset this environment variable because otherwise when we import `_import_snopt_from_path`
# the snopt module gets automatically imported, thus failing the import test below
os.environ.pop("PYOPTSPARSE_IMPORT_SNOPT_FROM", None)

# First party modules
from pyoptsparse.pySNOPT.pySNOPT import _import_snopt_from_path  # noqa: E402


class TestImportSnoptFromPath(unittest.TestCase):
    def test_nonexistent_path(self):
        with self.assertWarns(ImportWarning):
            self.assertIsNone(_import_snopt_from_path("/a/nonexistent/path"))

    def test_sys_path_unchanged(self):
        path = tuple(sys.path)
        _import_snopt_from_path("/some/path")
        self.assertEqual(tuple(sys.path), path)
