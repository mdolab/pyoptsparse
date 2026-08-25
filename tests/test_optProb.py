"""Test for optProb"""

# Standard Python modules
import unittest

# External modules
from baseclasses.testing.assertions import assert_dict_allclose, assert_dict_not_allclose, assert_not_allclose
import numpy as np
from numpy.testing import assert_allclose

try:
    # External modules
    from mpi4py import MPI

    comm = MPI.COMM_WORLD
except ImportError:
    comm = None

# First party modules
from pyoptsparse import OPT, Optimization
from pyoptsparse.pyOpt_utils import INFINITY, convertToCSR, convertToDense
from pyoptsparse.testing.pyOpt_testing import assert_optProb_size


class TestOptProb(unittest.TestCase):
    N_PROCS = 2
    tol = 1e-12

    def setUp(self):
        self.rng = np.random.default_rng(12345)

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
            lower = self.rng.uniform(-5, 2, n)
            upper = self.rng.uniform(5, 20, n)
            x0 = self.rng.uniform(lower, upper)
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
            lower = self.rng.uniform(-5, 2, nc)
            upper = self.rng.uniform(5, 6, nc)
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
            xScale=[self.rng.random(i) for i in nDV],
            objScale=[0.3],
            conScale=[self.rng.random(i) for i in nCon],
            offset=[self.rng.random(i) * np.arange(i) for i in nDV],
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
            """Helper function since some functions have optional arguments that are needed"""
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

    def test_parallel_add(self):
        """Check that when different procs add different variables/constraints/objectives, they are all collected
        properly when finalized.
        """
        if comm is None:
            raise unittest.SkipTest("mpi4py not available, skipping test.")
        if comm.size != self.N_PROCS:
            raise unittest.SkipTest("Not running with %d MPI procs, skipping test." % self.N_PROCS)
        objNames = [f"parallel-obj{i}" for i in range(comm.size)]
        dvNames = [f"parallel-x{i}" for i in range(comm.size)]
        conNames = [f"parallel-con{i}" for i in range(comm.size)]

        self.setup_optProb(nObj=1, nDV=[4, 8], nCon=[2, 3], xScale=[1.0, 1.0], conScale=[1.0, 1.0], offset=[0, 0])
        self.optProb.addObj(objNames[comm.rank])
        self.optProb.addVar(dvNames[comm.rank])
        self.optProb.addCon(conNames[comm.rank])
        self.optProb.finalize()

        # Get variables/constraints/objectives from each proc and check they are all the same
        allObjNames = comm.allgather(list(self.optProb.objectives.keys()))
        self.assertEqual(allObjNames[0], allObjNames[1])
        allDVNames = comm.allgather(list(self.optProb.variables.keys()))
        self.assertEqual(allDVNames[0], allDVNames[1])
        allConNames = comm.allgather(list(self.optProb.constraints.keys()))
        self.assertEqual(allConNames[0], allConNames[1])


class TestScaling(unittest.TestCase):
    def setUp(self):
        # Distinct, non-trivial per-element scales and offsets so that any
        # mixed-up indexing or row/column confusion shows up.
        self.xScale = {"x": [2.0, 0.5, 4.0], "y": [10.0, 0.1]}
        self.xOffset = {"x": [1.0, -2.0, 0.5], "y": [0.0, 3.0]}
        self.objScale = 3.0
        self.conScaleVals = {"c1": [5.0, 0.2], "c2": [7.0]}

        def objfunc(xdict):
            # Never actually called in these tests, but required by the API.
            return {"obj": 0.0, "c1": np.zeros(2), "c2": np.zeros(1)}, False

        optProb = Optimization("scaling-test", objfunc)
        optProb.addVarGroup("x", 3, lower=-10, upper=10, scale=self.xScale["x"], offset=self.xOffset["x"])
        optProb.addVarGroup("y", 2, lower=-10, upper=10, scale=self.xScale["y"], offset=self.xOffset["y"])
        optProb.addObj("obj", scale=self.objScale)
        optProb.addConGroup("c1", 2, lower=-1, upper=1, scale=self.conScaleVals["c1"])
        optProb.addConGroup("c2", 1, lower=-1, upper=1, scale=self.conScaleVals["c2"])
        optProb.finalize()

        self.optProb = optProb
        self.ndvs = optProb.ndvs
        self.nCon = optProb.nCon
        # invXScale = 1/scale, in DV order (x then y)
        self.invXScale = optProb.invXScale
        # conScale in natural (un-reordered) order: c1, c2
        self.conScale = optProb.conScale

    def test_finalize_scale(self):
        """Check that finalize() assembles per-DV/constraint scale and offset arrays
        in the correct flattened order.
        """
        assert_allclose(self.invXScale, 1.0 / np.array([2.0, 0.5, 4.0, 10.0, 0.1]))
        assert_allclose(self.conScale, [5.0, 0.2, 7.0])
        assert_allclose(self.optProb.xOffset, [1.0, -2.0, 0.5, 0.0, 3.0])

    def test_mapX_roundtrip(self):
        """Check that _mapXtoOpt matches the documented formula and that _mapXtoUser
        exactly inverts it.
        """
        rng = np.random.default_rng(0)
        x_user = rng.uniform(-5, 5, self.ndvs)
        x_opt = self.optProb._mapXtoOpt(x_user)
        # x_opt = (x_user - offset) / invXScale
        assert_allclose(x_opt, (x_user - self.optProb.xOffset) / self.invXScale)
        # round trip
        assert_allclose(self.optProb._mapXtoUser(x_opt), x_user)

    def test_mapObjGrad(self):
        """Check the objective gradient mapping g_opt = g_user * s_f * invXScale, and that
        the mapping does not mutate its input array in place.
        """
        # Objective gradient mapping: g_opt = g_user * s_f * invXScale (column/chain-rule scaling).
        rng = np.random.default_rng(1)
        gobj = rng.uniform(-3, 3, (self.optProb.nObj, self.ndvs))
        gobj_orig = gobj.copy()
        gobj_opt = self.optProb._mapObjGradtoOpt(gobj)
        assert_allclose(gobj_opt, gobj * self.objScale * self.invXScale)
        # the method must not mutate its input
        assert_allclose(gobj, gobj_orig)

    def test_mapConJac_roundtrip(self):
        """Check the in-place constraint Jacobian scaling J_opt = diag(conScale) . J . diag(invXScale),
        and that _mapConJactoUser inverts it back to the original values.
        """
        # Build an arbitrary dense Jacobian of the right shape and convert to CSR.
        rng = np.random.default_rng(2)
        dense = rng.uniform(-2, 2, (self.nCon, self.ndvs))
        jac = convertToCSR(dense)

        # _mapConJactoOpt works in place: J_opt = diag(conScale) . J . diag(invXScale)
        self.optProb._mapConJactoOpt(jac)
        expected = np.diag(self.conScale) @ dense @ np.diag(self.invXScale)
        assert_allclose(convertToDense(jac), expected)

        # _mapConJactoUser must invert it back to the original.
        self.optProb._mapConJactoUser(jac)
        assert_allclose(convertToDense(jac), dense)

    def test_mapObj_roundtrip(self):
        """Check the scalar objective value mapping and its round trip through
        _mapObjtoOpt / _mapObjtoUser.
        """
        f_user = 2.5
        f_opt = self.optProb._mapObjtoOpt(f_user)
        assert_allclose(f_opt, f_user * self.objScale)
        assert_allclose(self.optProb._mapObjtoUser(f_opt), f_user)

    def test_mapCon_roundtrip(self):
        """Check the constraint value mapping and its round trip through
        _mapContoOpt / _mapContoUser.
        """
        c_user = [1.0, -2.0, 3.0]
        c_opt = self.optProb._mapContoOpt(c_user)
        assert_allclose(c_opt, c_user * self.conScale)
        assert_allclose(self.optProb._mapContoUser(c_opt), c_user)

    def test_scale_and_offset(self):
        """A DV group with both a non-unit scale and a non-zero offset is the
        classic place to get the order of operations wrong.
        """

        def objfunc(xdict):
            return {"obj": 0.0}, False

        optProb = Optimization("edge", objfunc)
        optProb.addVarGroup("x", 2, lower=-10, upper=10, scale=4.0, offset=3.0)
        optProb.addObj("obj")
        optProb.finalize()

        x_user = [3.0, 7.0]  # note x_user[0] == offset
        x_opt = optProb._mapXtoOpt(x_user)
        # (x - 3) * 4
        assert_allclose(x_opt, [0.0, 16.0])
        assert_allclose(optProb._mapXtoUser(x_opt), x_user)

    def test_infinite_bounds_not_scaled(self):
        """INFINITY bounds must remain unbounded; scale/offset must not turn
        them into finite numbers in the assembled bounds.
        """

        def objfunc(xdict):
            return {"obj": 0.0}, False

        optProb = Optimization("inf", objfunc)
        optProb.addVarGroup("x", 1, lower=None, upper=None, scale=10.0, offset=5.0)
        optProb.addObj("obj")
        optProb.finalize()

        var = optProb.variables["x"][0]
        # Variable stores scaled bounds; unbounded sides stay at exactly +/- INFINITY.
        self.assertEqual(var.lower, -INFINITY)
        self.assertEqual(var.upper, INFINITY)


if __name__ == "__main__":
    unittest.main()
