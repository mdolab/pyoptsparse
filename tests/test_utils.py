"""
Unit tests for the sparse-matrix utilities in pyOpt_utils.

These functions (format conversions, index maps, row/column scaling, row
extraction) are used throughout the mapping/scaling layer but were previously
only exercised indirectly through full optimizations. Here we test them
directly on small hand-built matrices where the correct answer is obvious.
"""

# Standard Python modules
import unittest

# External modules
import numpy as np
from numpy.testing import assert_allclose, assert_array_equal
from scipy import sparse

# First party modules
from pyoptsparse.pyOpt_utils import (
    INFINITY,
    _broadcast_to_array,
    convertToCOO,
    convertToCSC,
    convertToCSR,
    convertToDense,
    extractRows,
    mapToCSC,
    mapToCSR,
    scaleColumns,
    scaleRows,
)


class TestSparseConversions(unittest.TestCase):
    def setUp(self):
        # A small, asymmetric, genuinely sparse reference matrix.
        # | 1 0 2 |
        # | 0 3 0 |
        # | 4 0 5 |
        self.dense = np.array(
            [
                [1.0, 0.0, 2.0],
                [0.0, 3.0, 0.0],
                [4.0, 0.0, 5.0],
            ]
        )
        # COO representation, intentionally NOT in row-major order to make sure
        # the conversion routines sort correctly.
        self.coo = {
            "coo": [
                np.array([2, 0, 1, 0, 2]),
                np.array([2, 0, 1, 2, 0]),
                np.array([5.0, 1.0, 3.0, 2.0, 4.0]),
            ],
            "shape": [3, 3],
        }

    def test_coo_to_dense(self):
        assert_allclose(convertToDense(self.coo), self.dense)

    def test_roundtrip_through_all_formats(self):
        # COO -> CSR -> CSC -> COO -> dense should reproduce the original.
        csr = convertToCSR(self.coo)
        csc = convertToCSC(csr)
        coo2 = convertToCOO(csc)
        assert_allclose(convertToDense(csr), self.dense)
        assert_allclose(convertToDense(csc), self.dense)
        assert_allclose(convertToDense(coo2), self.dense)

    def test_convertToCOO_passthrough(self):
        # Already-COO input should be returned unchanged.
        self.assertIs(convertToCOO(self.coo), self.coo)

    def test_convertToCSR_idempotent(self):
        csr = convertToCSR(self.coo)
        self.assertIs(convertToCSR(csr), csr)

    def test_dense_array_input(self):
        # A plain numpy array should be accepted and converted correctly.
        assert_allclose(convertToDense(convertToCSR(self.dense)), self.dense)

    def test_unknown_format_raises(self):
        # A ragged nested list cannot be coerced into a dense array, so it
        # should fall through to the explicit ValueError.
        with self.assertRaises(ValueError):
            convertToCOO([[1.0, 2.0], [3.0]])


class TestIndexMaps(unittest.TestCase):
    """mapToCSR/mapToCSC return index arrays into the original data; the subtle
    part is that the permutation must reproduce the matrix exactly."""

    def setUp(self):
        self.dense = np.array(
            [
                [1.0, 0.0, 2.0],
                [0.0, 3.0, 0.0],
                [4.0, 0.0, 5.0],
            ]
        )
        self.coo = {
            "coo": [
                np.array([2, 0, 1, 0, 2]),
                np.array([2, 0, 1, 2, 0]),
                np.array([5.0, 1.0, 3.0, 2.0, 4.0]),
            ],
            "shape": [3, 3],
        }

    def test_mapToCSR(self):
        row_p, col_idx, idx_data = mapToCSR(self.coo)
        data = np.asarray(self.coo["coo"][2])
        csr = {"csr": [row_p, col_idx, data[idx_data]], "shape": [3, 3]}
        assert_allclose(convertToDense(csr), self.dense)
        # last entry of the row pointer is the nnz
        self.assertEqual(row_p[-1], 5)

    def test_mapToCSC(self):
        row_idx, col_p, idx_data = mapToCSC(self.coo)
        data = np.asarray(self.coo["coo"][2])
        csc = {"csc": [col_p, row_idx, data[idx_data]], "shape": [3, 3]}
        assert_allclose(convertToDense(csc), self.dense)
        self.assertEqual(col_p[-1], 5)


class TestRowColScaling(unittest.TestCase):
    def setUp(self):
        self.dense = np.array(
            [
                [1.0, 0.0, 2.0],
                [0.0, 3.0, 0.0],
                [4.0, 0.0, 5.0],
            ]
        )

    def test_scaleRows(self):
        csr = convertToCSR(self.dense)
        factor = np.array([10.0, 100.0, 1000.0])
        scaleRows(csr, factor)
        assert_allclose(convertToDense(csr), np.diag(factor) @ self.dense)

    def test_scaleColumns(self):
        csr = convertToCSR(self.dense)
        factor = np.array([2.0, 3.0, 4.0])
        scaleColumns(csr, factor)
        assert_allclose(convertToDense(csr), self.dense @ np.diag(factor))

    def test_scale_wrong_length_raises(self):
        csr = convertToCSR(self.dense)
        with self.assertRaises(ValueError):
            scaleRows(csr, np.array([1.0, 2.0]))
        with self.assertRaises(ValueError):
            scaleColumns(csr, np.array([1.0, 2.0]))

    def test_scale_requires_csr(self):
        coo = convertToCOO(self.dense)
        with self.assertRaises(ValueError):
            scaleRows(coo, np.array([1.0, 1.0, 1.0]))


class TestExtractRows(unittest.TestCase):
    def test_extractRows(self):
        dense = np.array(
            [
                [1.0, 0.0, 2.0],
                [0.0, 3.0, 0.0],
                [4.0, 0.0, 5.0],
            ]
        )
        csr = convertToCSR(dense)
        sub = extractRows(csr, [0, 2])
        assert_allclose(convertToDense(sub), dense[[0, 2], :])
        self.assertEqual(sub["shape"], [2, 3])


class TestScipySparseWarning(unittest.TestCase):
    def test_scipy_input_warns(self):
        spmat = sparse.csr_matrix(np.array([[1.0, 0.0], [0.0, 2.0]]))
        with self.assertWarns(UserWarning):
            coo = convertToCOO(spmat)
        assert_allclose(convertToDense(coo), np.array([[1.0, 0.0], [0.0, 2.0]]))


class TestBroadcastToArray(unittest.TestCase):
    def test_scalar_broadcast(self):
        out = _broadcast_to_array("scale", 2.0, 4)
        assert_array_equal(out, np.full(4, 2.0))

    def test_array_passthrough(self):
        out = _broadcast_to_array("scale", [1.0, 2.0, 3.0], 3)
        assert_array_equal(out, np.array([1.0, 2.0, 3.0]))

    def test_wrong_length_raises(self):
        with self.assertRaises(ValueError):
            _broadcast_to_array("scale", [1.0, 2.0], 3)

    def test_none_disallowed_by_default(self):
        with self.assertRaises(ValueError):
            _broadcast_to_array("lower", None, 3)

    def test_none_allowed_when_requested(self):
        out = _broadcast_to_array("lower", None, 3, allow_none=True)
        self.assertEqual(len(out), 3)
        self.assertTrue(all(v is None for v in out))


class TestConstants(unittest.TestCase):
    def test_infinity_value(self):
        self.assertEqual(INFINITY, 1e20)


if __name__ == "__main__":
    unittest.main()
