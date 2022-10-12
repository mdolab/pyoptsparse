"""Test for optProb"""

# Standard Python modules
import unittest

# External modules
from baseclasses.testing.assertions import assert_dict_allclose, assert_dict_not_allclose, assert_not_allclose
import numpy as np
from numpy.testing import assert_allclose

# First party modules
from pyoptsparse import OPT, Optimization

# Local modules
from testing_utils import assert_optProb_size


class TestOptProb(unittest.TestCase):
    tol = 1e-12

    def objfunc(self, xdict):
        """
        This is a simple quadratic test function with linear constraints.
        The actual problem doesn't really matter, since we are not testing optimization,
        but just optProb. However, we need to initialize and run an optimization
        in order to have optimizer-specific fields in optProb populated, such as
        jacIndices.

        This problem is probably not feasible, but that's okay.
        """
        funcs = {}
        funcs["obj_0"] = 0
        for x in xdict.keys():
            funcs["obj_0"] += np.sum(np.power(xdict[x], 2))
        for iCon, nc in enumerate(self.nCon):
            conName = f"con_{iCon}"
            funcs[conName] = np.zeros(nc)
            for x in xdict.keys():
                for j in range(nc):
                    funcs[conName][j] = (iCon + 1) * np.sum(xdict[x])
        return funcs, False

    def setup_optProb(self, nObj=1, nDV=[4], nCon=[2], xScale=[1.0], objScale=[1.0], conScale=[1.0], offset=[0.0]):
        """
        This function sets up a general optimization problem, with arbitrary
        DVs, constraints and objectives.
        Arbitrary scaling for the various parameters can also be specified.
        """
        self.nObj = nObj
        self.nDV = nDV
        self.nCon = nCon
        self.xScale = xScale
        self.objScale = objScale
        self.conScale = conScale
        self.offset = offset

        # Optimization Object
        self.optProb = Optimization("Configurable Test Problem", self.objfunc)
        self.x0 = {}
        # Design Variables
        for iDV in range(len(nDV)):
            n = nDV[iDV]
            lower = np.random.uniform(-5, 2, n)
            upper = np.random.uniform(5, 20, n)
            x0 = np.random.uniform(lower, upper)
            dvName = f"x{iDV}"
            self.x0[dvName] = x0
            self.optProb.addVarGroup(
                dvName,
                n,
                lower=lower,
                upper=upper,
                value=x0,
                scale=xScale[iDV],
                offset=offset[iDV],
            )

        # Constraints
        for iCon in range(len(nCon)):
            nc = nCon[iCon]
            lower = np.random.uniform(-5, 2, nc)
            upper = np.random.uniform(5, 6, nc)
            self.optProb.addConGroup(
                f"con_{iCon}",
                nc,
                lower=lower,
                upper=upper,
                scale=conScale[iCon],
            )

        # Objective
        for iObj in range(nObj):
            self.optProb.addObj(f"obj_{iObj}", scale=objScale[iObj])

        # Finalize
        self.optProb.printSparsity()
        # run optimization
        # we don't care about outputs, but this performs optimizer-specific re-ordering
        # of constraints so we need this to test mappings
        opt = OPT("slsqp", options={"IFILE": "optProb_SLSQP.out"})
        opt(self.optProb, "FD")

    def test_setDV_getDV(self):
        """
        We just test that setDV and getDV work, even with scaling
        """
        self.setup_optProb(
            nObj=1,
            nDV=[4, 8],
            nCon=[2, 3],
            xScale=[4, 1],
            objScale=[0.3],
            conScale=[0.1, 8],
            offset=[3, 7],
        )
        # test getDV first
        x0 = self.optProb.getDVs()
        assert_dict_allclose(x0, self.x0)
        # now set, get, and compare
        newDV = {"x0": np.arange(4), "x1": np.arange(8)}
        self.optProb.setDVs(newDV)
        outDV = self.optProb.getDVs()
        assert_dict_allclose(newDV, outDV)

    def test_setDV_VarGroup(self):
        """
        Test that setDV works with a subset of VarGroups
        """
        self.setup_optProb(
            nObj=1,
            nDV=[4, 8],
            nCon=[2, 3],
            xScale=[4, 1],
            objScale=[0.3],
            conScale=[0.1, 8],
            offset=[3, 7],
        )
        oldDV = self.optProb.getDVs()
        # set values for only one VarGroup
        newDV = {"x0": np.arange(4)}
        self.optProb.setDVs(newDV)
        outDV = self.optProb.getDVs()
        # check x0 changed
        assert_allclose(newDV["x0"], outDV["x0"])
        # check x1 is the same
        assert_allclose(oldDV["x1"], outDV["x1"])

    def test_mappings(self):
        """
        This test checks the various mapping and process helper functions
        in pyOpt_optimization. In this function we just set up an optimization problem,
        and the actual test is done in `map_check_value`.
        """
        nDV = [4, 8, 1]
        nCon = [2, 3, 1, 1]
        self.setup_optProb(
            nObj=1,
            nDV=nDV,
            nCon=nCon,
            xScale=[np.random.rand(i) for i in nDV],
            objScale=[0.3],
            conScale=[np.random.rand(i) for i in nCon],
            offset=[np.random.rand(i) * np.arange(i) for i in nDV],
        )

        # first test X
        x = self.optProb.getDVs()
        self.map_check_value("X", x)

        # next we check the objective
        funcs, _ = self.objfunc(x)
        obj_funcs = {}
        for key in funcs.keys():
            if "obj" in key:
                obj_funcs[key] = funcs[key]
        self.map_check_value("Obj", obj_funcs)

        # lastly we check the constraints
        funcs, _ = self.objfunc(x)
        con_funcs = {}
        for key in funcs.keys():
            if "con" in key:
                con_funcs[key] = funcs[key]
        self.map_check_value("Con", con_funcs)

    def map_check_value(self, key, val):
        """
        This function checks all the mapping and process functions
        in both directions, for a given key = {'X', 'Con', 'Obj'}
        and val in dictionary format.
        """
        # dictionary of function handles to test
        map_funcs = {
            "X": [self.optProb._mapXtoOpt, self.optProb._mapXtoUser],
            "X_Dict": [self.optProb._mapXtoOpt_Dict, self.optProb._mapXtoUser_Dict],
            "Con": [self.optProb._mapContoOpt, self.optProb._mapContoUser],
            "Con_Dict": [self.optProb._mapContoOpt_Dict, self.optProb._mapContoUser_Dict],
            "Obj": [self.optProb._mapObjtoOpt, self.optProb._mapObjtoUser],
            "Obj_Dict": [self.optProb._mapObjtoOpt_Dict, self.optProb._mapObjtoUser_Dict],
        }
        process_funcs = {
            "X": {"vec": self.optProb.processXtoVec, "dict": self.optProb.processXtoDict},
            "Con": {"vec": self.optProb.processContoVec, "dict": self.optProb.processContoDict},
            "Obj": {"vec": self.optProb.processObjtoVec, "dict": self.optProb.processObjtoDict},
        }

        def processValue(key, val, output):
            """helper function since some functions have optional arguments that are needed"""
            if key == "Con":
                return process_funcs[key][output](val, scaled=False, natural=True)
            elif key == "Obj":
                return process_funcs[key][output](val, scaled=False)
            else:
                return process_funcs[key][output](val)

        # test dict to vec mappings
        vec = processValue(key, val, "vec")
        dictionary = processValue(key, vec, "dict")
        assert_dict_allclose(val, dictionary)

        # test mappings using dictionaries
        val_opt = map_funcs[key + "_Dict"][0](val)
        val_user = map_funcs[key + "_Dict"][1](val_opt)
        assert_dict_allclose(val_user, val)
        assert_dict_not_allclose(val_user, val_opt)

        # test mappings using vectors
        val = processValue(key, val, "vec")
        val_opt = map_funcs[key][0](val)
        val_user = map_funcs[key][1](val_opt)
        assert_allclose(val_user, val, atol=self.tol, rtol=self.tol)
        assert_not_allclose(val_user, val_opt)

        # check that the scaling was actually done correctly
        # we only check this for the array version because
        # it's much simpler
        if key == "X":
            scale = np.hstack(self.xScale)
            offset = np.hstack(self.offset)
            assert_allclose(val_opt, (val_user - offset) * scale)
        else:
            if key == "Obj":
                scale = np.hstack(self.objScale)
            else:
                scale = np.hstack(self.conScale)
            assert_allclose(val_opt, val_user * scale)

    def test_finalize(self):
        """
        Check that multiple finalize calls don't mess up the optProb
        """
        self.setup_optProb(nObj=1, nDV=[4, 8], nCon=[2, 3], xScale=[1.0, 1.0], conScale=[1.0, 1.0], offset=[0, 0])
        assert_optProb_size(self.optProb, 1, 12, 5)
        self.optProb.addObj("obj2")
        assert_optProb_size(self.optProb, 2, 12, 5)
        self.optProb.addVar("DV2")
        assert_optProb_size(self.optProb, 2, 13, 5)
        self.optProb.addCon("CON2")
        assert_optProb_size(self.optProb, 2, 13, 6)


if __name__ == "__main__":
    unittest.main()
