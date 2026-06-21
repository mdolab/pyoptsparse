"""
Unit tests for the sparse-matrix utilities in pyOpt_utils.
"""

# Standard Python modules
import unittest

# External modules
import numpy as np
from numpy.testing import assert_allclose, assert_array_equal

# First party modules
from pyoptsparse.pyOpt_utils import (
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

# All three sparse fixtures represent the same 3x3 matrix:
#   [[1, 0, 2],
#    [0, 3, 0],
#    [4, 0, 5]]
#
# _COO is intentionally NOT in row-major order to verify that conversion
# routines handle unsorted input correctly.
#
# _CSR is the expected output of convertToCSR(_COO): elements within each
# row appear in COO-arrival order, not sorted by column.
#   row 0: (col 0, 1.0), (col 2, 2.0)
#   row 1: (col 1, 3.0)
#   row 2: (col 2, 5.0), (col 0, 4.0)  <- arrival order from _COO
#
# _CSC is the expected output of convertToCSC(_CSR): elements within each
# column appear in row-scan order from the CSR pass.
#   col 0: (row 0, 1.0), (row 2, 4.0)
#   col 1: (row 1, 3.0)
#   col 2: (row 0, 2.0), (row 2, 5.0)
_DENSE = np.array(
    [
        [1.0, 0.0, 2.0],
        [0.0, 3.0, 0.0],
        [4.0, 0.0, 5.0],
    ]
)
_COO = {
    "coo": [
        [2, 0, 1, 0, 2],
        [2, 0, 1, 2, 0],
        [5.0, 1.0, 3.0, 2.0, 4.0],
    ],
    "shape": [3, 3],
}
_CSR = {
    "csr": [
        [0, 2, 3, 5],
        [0, 2, 1, 2, 0],
        [1.0, 2.0, 3.0, 5.0, 4.0],
    ],
    "shape": [3, 3],
}
_CSC = {
    "csc": [
        [0, 2, 3, 5],
        [0, 2, 1, 0, 2],
        [1.0, 4.0, 3.0, 2.0, 5.0],
    ],
    "shape": [3, 3],
}


class TestConvertToDense(unittest.TestCase):
    def test_from_coo(self):
        assert_allclose(convertToDense(_COO), _DENSE)

    def test_from_csr(self):
        assert_allclose(convertToDense(_CSR), _DENSE)

    def test_from_csc(self):
        assert_allclose(convertToDense(_CSC), _DENSE)

    def test_from_dense_array(self):
        assert_allclose(convertToDense(_DENSE), _DENSE)


class TestConvertToCOO(unittest.TestCase):
    def test_from_csr(self):
        coo = convertToCOO(_CSR)
        rows, cols, data = coo["coo"]
        assert_array_equal(rows, [0, 0, 1, 2, 2])
        assert_array_equal(cols, [0, 2, 1, 2, 0])
        assert_allclose(data, [1.0, 2.0, 3.0, 5.0, 4.0])
        self.assertEqual(coo["shape"], [3, 3])

    def test_from_csc(self):
        # CSC scans the column pointer, so COO entries emerge in column order, not _COO arrival order.
        # col 0: (row 0, 1.0),(row 2, 4.0); col 1: (row 1, 3.0); col 2: (row 0, 2.0),(row 2, 5.0)
        coo = convertToCOO(_CSC)
        rows, cols, data = coo["coo"]
        assert_array_equal(rows, [0, 2, 1, 0, 2])
        assert_array_equal(cols, [0, 0, 1, 2, 2])
        assert_allclose(data, [1.0, 4.0, 3.0, 2.0, 5.0])
        self.assertEqual(coo["shape"], [3, 3])

    def test_from_dense_array(self):
        coo = convertToCOO(_DENSE)
        rows, cols, data = coo["coo"]
        reconstructed = np.zeros(_DENSE.shape)
        for r, c, v in zip(rows, cols, data, strict=True):
            reconstructed[r, c] = v
        assert_allclose(reconstructed, _DENSE)

    def test_idempotent(self):
        self.assertIs(convertToCOO(_COO), _COO)


class TestConvertToCSR(unittest.TestCase):
    def test_from_coo(self):
        # COO arrives unordered; elements land in row buckets in COO-arrival order.
        # row 0: (col 0, 1.0) then (col 2, 2.0); row 2: (col 2, 5.0) then (col 0, 4.0)
        csr = convertToCSR(_COO)
        rowp, col_idx, data = csr["csr"]
        assert_array_equal(rowp, [0, 2, 3, 5])
        assert_array_equal(col_idx, [0, 2, 1, 2, 0])
        assert_allclose(data, [1.0, 2.0, 3.0, 5.0, 4.0])
        self.assertEqual(csr["shape"], [3, 3])

    def test_from_csc(self):
        # _CSC expands to COO in column-scan order; that COO then feeds the CSR builder.
        # row 0: (col 0, 1.0) then (col 2, 2.0); row 2: (col 0, 4.0) then (col 2, 5.0)
        csr = convertToCSR(_CSC)
        rowp, col_idx, data = csr["csr"]
        assert_array_equal(rowp, [0, 2, 3, 5])
        assert_array_equal(col_idx, [0, 2, 1, 0, 2])
        assert_allclose(data, [1.0, 2.0, 3.0, 4.0, 5.0])
        self.assertEqual(csr["shape"], [3, 3])

    def test_from_dense_array(self):
        csr = convertToCSR(_DENSE)
        self.assertIn("csr", csr)
        assert_allclose(convertToDense(csr), _DENSE)

    def test_idempotent(self):
        self.assertIs(convertToCSR(_CSR), _CSR)


class TestConvertToCSC(unittest.TestCase):
    def test_from_coo(self):
        # Converts COO -> CSR -> CSC. Column-scan order from _CSR:
        # col 0: (row 0, 1.0),(row 2, 4.0); col 1: (row 1, 3.0); col 2: (row 0, 2.0),(row 2, 5.0)
        csc = convertToCSC(_COO)
        colp, row_idx, data = csc["csc"]
        assert_array_equal(colp, [0, 2, 3, 5])
        assert_array_equal(row_idx, [0, 2, 1, 0, 2])
        assert_allclose(data, [1.0, 4.0, 3.0, 2.0, 5.0])
        self.assertEqual(csc["shape"], [3, 3])

    def test_from_csr(self):
        # _CSR -> CSC: same column-scan order, same result as test_from_coo.
        csc = convertToCSC(_CSR)
        colp, row_idx, data = csc["csc"]
        assert_array_equal(colp, [0, 2, 3, 5])
        assert_array_equal(row_idx, [0, 2, 1, 0, 2])
        assert_allclose(data, [1.0, 4.0, 3.0, 2.0, 5.0])
        self.assertEqual(csc["shape"], [3, 3])

    def test_from_dense_array(self):
        csc = convertToCSC(_DENSE)
        self.assertIn("csc", csc)
        assert_allclose(convertToDense(csc), _DENSE)

    def test_idempotent(self):
        self.assertIs(convertToCSC(_CSC), _CSC)


class TestIndexMaps(unittest.TestCase):
    """mapToCSR/mapToCSC return index arrays into the original data; the subtle
    part is that the permutation must reproduce the matrix exactly.
    """

    def test_mapToCSR(self):
        row_p, col_idx, idx_data = mapToCSR(_COO)
        data = np.asarray(_COO["coo"][2])
        csr = {"csr": [row_p, col_idx, data[idx_data]], "shape": [3, 3]}
        assert_allclose(convertToDense(csr), _DENSE)
        # last entry of the row pointer is the nnz
        self.assertEqual(row_p[-1], 5)

    def test_mapToCSC(self):
        row_idx, col_p, idx_data = mapToCSC(_COO)
        data = np.asarray(_COO["coo"][2])
        csc = {"csc": [col_p, row_idx, data[idx_data]], "shape": [3, 3]}
        assert_allclose(convertToDense(csc), _DENSE)
        self.assertEqual(col_p[-1], 5)


class TestRowColScaling(unittest.TestCase):
    def test_scaleRows(self):
        csr = convertToCSR(_DENSE)
        factor = np.array([10.0, 100.0, 1000.0])
        scaleRows(csr, factor)
        assert_allclose(convertToDense(csr), np.diag(factor) @ _DENSE)

    def test_scaleColumns(self):
        csr = convertToCSR(_DENSE)
        factor = np.array([2.0, 3.0, 4.0])
        scaleColumns(csr, factor)
        assert_allclose(convertToDense(csr), _DENSE @ np.diag(factor))

    def test_scale_wrong_length_raises(self):
        csr = convertToCSR(_DENSE)
        with self.assertRaises(ValueError):
            scaleRows(csr, np.array([1.0, 2.0]))
        with self.assertRaises(ValueError):
            scaleColumns(csr, np.array([1.0, 2.0]))


class TestExtractRows(unittest.TestCase):
    def test_extractRows(self):
        csr = convertToCSR(_DENSE)
        sub = extractRows(csr, [0, 2])
        assert_allclose(convertToDense(sub), _DENSE[[0, 2], :])
        self.assertEqual(sub["shape"], [2, 3])


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


if __name__ == "__main__":
    unittest.main()
