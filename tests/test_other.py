# Standard Python modules
import sys
import unittest

# First party modules
from pyoptsparse.pySNOPT.pySNOPT import _import_snopt_from_path


class TestImportSnoptFromPath(unittest.TestCase):
    def test_nonexistent_path(self):
        assert _import_snopt_from_path("/a/nonexistent/path") is None

    def test_sys_path_unchanged(self):
        path = tuple(sys.path)
        _import_snopt_from_path("/some/path")
        assert tuple(sys.path) == path
