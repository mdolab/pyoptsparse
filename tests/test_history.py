"""Test History.getValues with major and minor iterations (issues #250 and #269)"""

# Standard Python modules
import unittest

# External modules
from sqlitedict import SqliteDict

# First party modules
from pyoptsparse import OPT, History, Optimization


def objfunc(xdict):
    x = xdict["xvars"]
    funcs = {"obj": (x[0] - 1) ** 2 + (x[1] - 2) ** 2, "con": [x[0] + x[1]]}
    return funcs, False


class TestHistoryMajorIterations(unittest.TestCase):
    def setUp(self):
        self.histFileName = f"history_{self._testMethodName}.hst"

        optProb = Optimization("quadratic", objfunc)
        optProb.addVarGroup("xvars", 2, lower=-5, upper=5, value=0.0)
        optProb.addConGroup("con", 1, upper=2.0)
        optProb.addObj("obj")
        options = {"SwarmSize": 4, "maxOuterIter": 2, "maxInnerIter": 2, "minInnerIter": 1, "seed": 1, "fileout": 0}
        OPT("ALPSO", options=options)(optProb, storeHistory=self.histFileName)

        # Make every second function evaluation a minor iteration, and store values at the
        # major iterations only, as SNOPT does
        with SqliteDict(self.histFileName, autocommit=True) as db:
            callCounters = sorted((key for key in db.keys() if key.isdigit()), key=int)
            self.nMajor = 0
            for i, key in enumerate(callCounters):
                val = db[key]
                if i % 2 == 1:
                    val["isMajor"] = False
                else:
                    val.update({"isMajor": True, "nMajor": self.nMajor, "merit": float(i)})
                    self.nMajor += 1
                db[key] = val
        self.nCalls = len(callCounters)

    def test_value_stored_at_all_iterations(self):
        hist = History(self.histFileName)
        data = hist.getValues(names="isMajor", major=False)
        self.assertEqual(data["isMajor"].shape, (self.nCalls, 1))
        self.assertEqual(int(data["isMajor"].sum()), self.nMajor)

    def test_all_names_minor_iterations(self):
        hist = History(self.histFileName)
        data = hist.getValues(major=False)
        self.assertNotIn("merit", data)
        self.assertNotIn("nMajor", data)
        self.assertEqual(data["obj"].shape, (self.nCalls, 1))
        self.assertEqual(data["xvars"].shape, (self.nCalls, 2))

    def test_major_only_name_sets_major(self):
        hist = History(self.histFileName)
        data = hist.getValues(names=["obj", "merit"], major=False)
        self.assertEqual(data["merit"].shape, (self.nMajor, 1))
        self.assertEqual(data["obj"].shape, (self.nMajor, 1))

    def test_major_iterations(self):
        hist = History(self.histFileName)
        data = hist.getValues(major=True)
        self.assertEqual(data["merit"].shape, (self.nMajor, 1))
        self.assertEqual(data["obj"].shape, (self.nMajor, 1))


if __name__ == "__main__":
    unittest.main()
