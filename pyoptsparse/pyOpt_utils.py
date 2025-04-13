"""
pyOptSparse_utils holds a minimal set of sparse-matrix type routines for pyOptSparse.
This is designed to replace the SciPy sparse matrix formats, which have no way to enforce
a constant sparsity structure as required by the optimizers.
We use a very simple dictionary format to represent the three most common forms of sparse matrices::

    mat = {'coo':[row,  col,    data], 'shape':[nrow, ncols]} # A coo matrix
    mat = {'csr':[rowp, colind, data], 'shape':[nrow, ncols]} # A csr matrix
    mat = {'csc':[colp, rowind, data], 'shape':[nrow, ncols]} # A csc matrix
"""

# Standard Python modules
import importlib
import os
import sys
import types
from typing import Optional, Tuple, Union
import warnings

# External modules
import numpy as np
from numpy import ndarray
from scipy import sparse
from scipy.sparse import spmatrix

# Local modules
from .pyOpt_types import ArrayType

# Define index mnemonics
IROW = 0
ICOL = 1

IROWP = 0
ICOLIND = 1

ICOLP = 0
IROWIND = 1

IDATA = 2

# Constants
INFINITY = 1e20
EPS = np.finfo(np.float64).eps


def mapToCSR(mat: dict) -> Tuple[ndarray, ndarray, ndarray]:
    """
    Given a pyoptsparse matrix definition, return a tuple containing a
    map of the matrix to the CSR format.

    Parameters
    ----------
    mat : dict
       A sparse matrix representation.

    Returns
    -------
    tup : tuple of numpy arrays
        tup[0] : numpy array (size=num_rows+1)
            An array that holds the indices in col_idx and data at which each
            row begins.  The last index of contains the number of nonzero
            elements in the sparse array.
        tup[1] : numpy array (size=nnz)
            An array of the column indices of each element in data.
        tup[2] : numpy array (size=nnz)
            An indexing array which maps the elements in the data array
            to elements in the CSR data array.
    """
    if "csr" in mat:
        # First handle the trivial case CSR->CSR
        row_p = mat["csr"][IROW]
        col_idx = mat["csr"][ICOL]
        idx_data = np.s_[:]
        return row_p, col_idx, idx_data

    num_rows = mat["shape"][0]
    num_cols = mat["shape"][1]

    if "csc" in mat:
        # If given a CSC matrix, expand the column pointers so we
        # effectively have a COO representation.
        csc_colp = mat["csr"][ICOL]
        rows = mat["csc"][IROW]
        nnz = csc_colp[-1]

        # Allocate the COO maps
        cols = np.zeros(nnz, dtype="intc")

        # We already have a full representation of the columns.
        # We need to decompress the representation of the rows.
        for j in range(num_cols):
            cols[csc_colp[j] : csc_colp[j + 1]] = j

    elif "coo" in mat:
        rows = mat["coo"][IROW]
        cols = mat["coo"][ICOL]
        nnz = len(rows)

    # Allocate the row pointer array
    row_p = np.zeros(num_rows + 1, dtype="intc")

    # Get the sort order that puts data in row-major form
    idx_data = np.lexsort((cols, rows))

    # Apply the row-major indexing to the COO column and row indices
    col_idx = np.asarray(cols, dtype="intc")[idx_data]
    rows_rowmaj = np.asarray(rows, dtype="intc")[idx_data]

    # Now for i = 0 to num_rows-1, row_p[i] is the first occurrence
    # of i in rows_rowmaj
    row_p[:-1] = np.digitize(np.arange(num_rows), rows_rowmaj, right=True)

    # By convention store nnz in the last element of row_p
    row_p[-1] = nnz

    return row_p, col_idx, idx_data


def mapToCSC(mat: dict) -> Tuple[ndarray, ndarray, ndarray]:
    """
    Given a pyoptsparse matrix definition, return a tuple containing a
    map of the matrix to the CSC format.

    Parameters
    ----------
    mat : dict
       A sparse matrix representation.

    Returns
    -------
    tup : tuple of numpy arrays
        tup[0] : numpy array (size=nnz)
            An array that holds the row index of each element in the CSC
            representation of the data.
        tup[1] : numpy array (size=num_cols+1)
            An array that holds the indices in the CSC representation
            and data at which each column begins.  The last index of
            contains the number of nonzero elements in the sparse array.
        tup[2] : numpy array
            An indexing array which maps the elements in the data array
            to elements in the CSC data array.
    """
    if "csc" in mat:
        # First handle the trivial case CSR->CSR
        row_idx = mat["csc"][IROW]
        col_p = mat["csc"][ICOL]
        idx_data = np.s_[:]
        return row_idx, col_p, idx_data

    num_rows = mat["shape"][0]
    num_cols = mat["shape"][1]

    if "csr" in mat:
        # If given a CSR matrix, expand the row pointers so we
        # effectively have a COO representation.
        csr_rowp = mat["csr"][IROW]
        cols = mat["csr"][ICOL]
        nnz = csr_rowp[-1]

        # Allocate the COO maps
        rows = np.zeros(nnz, dtype="intc")

        # We already have a full representation of the columns.
        # We need to decompress the representation of the rows.
        for j in range(num_rows):
            rows[csr_rowp[j] : csr_rowp[j + 1]] = j

        # Now we have rows and cols, proceed as if we started with a COO matrix

    elif "coo" in mat:
        rows = mat["coo"][IROW]
        cols = mat["coo"][ICOL]
        nnz = len(rows)

    else:
        raise ValueError("Invalid matrix type")

    # Allocate the new column pointer
    col_p = np.zeros(num_cols + 1, dtype="intc")

    # Get the sort order that puts data in column-major form
    idx_data = np.lexsort((rows, cols))

    # Apply the column-major indexing to the COO column and row indices
    row_idx = np.asarray(rows, dtype="intc")[idx_data]
    cols_colmaj = np.asarray(cols, dtype="intc")[idx_data]

    # Now for i = 0 to num_cols-1, col_p[i] is the first occurrence
    # of i in cols_colmaj
    col_p[:-1] = np.digitize(np.arange(num_cols), cols_colmaj, right=True)

    # By convention store nnz in the last element of col_p
    col_p[-1] = nnz

    return row_idx, col_p, idx_data


def convertToCOO(mat: Union[dict, spmatrix, ndarray]):
    """
    Take a pyoptsparse sparse matrix definition of a COO, CSR or
    CSC matrix or numpy array or scipy sparse matrix and return
    the same matrix in COO format.

    Parameters
    ----------
    mat : dict or numpy array
       A sparse matrix representation or numpy array

    Returns
    -------
    newMat : dict
        A coo representation of the same matrix
    """

    if isinstance(mat, dict):
        if "coo" in mat:
            return mat
        if "csr" in mat:
            return _csr_to_coo(mat)
        elif "csc" in mat:
            return _csc_to_coo(mat)
    else:
        # Try to do it with a scipy sparse matrix:
        try:
            if sparse.issparse(mat):
                warnings.warn(
                    "Using scipy.sparse matrices with pyOptSparse is VERY STRONGLY discouraged. "
                    + "Please use the simplified pyOptSparse format which allows for "
                    + "fixed sparsity structure and explicit zeros in the matrix. "
                    + "There is no way to guarantee a fixed sparsity structure with scipy matrices "
                    + "which is what the underlying optimizers require. "
                    + "Using scipy.sparse matrices may cause unexpected errors.",
                    stacklevel=2,
                )

                mat = mat.tocoo()
                return {"coo": [mat.row, mat.col, mat.data], "shape": mat.shape}
        except Exception:
            pass

        # Now try to do it with a numpy matrix:
        try:
            return _denseToCOO(np.atleast_2d(np.array(mat)))
        except Exception as e:
            raise ValueError(
                "Unknown matrix format. "
                + "Must be a dense numpy array or a pyOptSparse sparse matrix format of COO, CSR or CSC. "
                + f"See documentation for correct format. Supplied Matrix is: {repr(mat)}"
            ) from e


def convertToCSR(mat: Union[dict, spmatrix, ndarray]) -> dict:
    """
    Take a pyoptsparse sparse matrix definition of a COO, CSR or
    CSC matrix or numpy array and return the same matrix in CSR format

    Parameters
    ----------
    mat : dict or numpy array
       A sparse matrix representation or numpy array

    Returns
    -------
    newMat : dict
        A coo representation of the same matrix
    """
    if isinstance(mat, dict) and "csr" in mat:
        return mat

    mat = convertToCOO(mat)
    n = mat["shape"][0]
    m = mat["shape"][1]
    rows = mat["coo"][IROW]
    cols = mat["coo"][ICOL]
    data = mat["coo"][IDATA]

    rowp = np.zeros(n + 1, dtype="intc")

    # Count up the number of times things are index
    for row in rows:
        rowp[row + 1] += 1

    # Set up the array as a pointer
    for i in range(1, n + 1):
        rowp[i] += rowp[i - 1]

    ncols = np.zeros(rowp[-1], dtype="intc")
    ndata = np.zeros(rowp[-1], dtype=type(data[0]))

    # Now, add all the values and the data
    for i in range(len(rows)):
        r = rows[i]
        ncols[rowp[r]] = cols[i]
        ndata[rowp[r]] = data[i]
        rowp[r] += 1

    # Readjust the pointer
    for i in range(n, 0, -1):
        rowp[i] = rowp[i - 1]
    rowp[0] = 0

    return {"csr": [rowp, ncols, ndata], "shape": [n, m]}


def convertToCSC(mat: Union[dict, spmatrix, ndarray]) -> dict:
    """
    Take a pyoptsparse sparse matrix definition of a COO, CSR or
    CSC matrix or numpy array and return the same matrix in CSR format

    Parameters
    ----------
    mat : dict or numpy array
       A sparse matrix representation or numpy array

    Returns
    -------
    newMat : dict
        A coo representation of the same matrix
    """
    if "csc" in mat:
        return mat

    mat = convertToCSR(mat)
    n = mat["shape"][0]
    m = mat["shape"][1]
    rowp = mat["csr"][IROWP]
    cols = mat["csr"][ICOLIND]
    data = mat["csr"][IDATA]

    # Allocate the new arrays
    colp = np.zeros(m + 1, "intc")
    rows = np.zeros(len(cols), "intc")

    # Count up the number of references to each column
    for col in cols:
        colp[col + 1] += 1

    # Set colp so that it is now a pointer
    for i in range(1, m):
        colp[i] += colp[i - 1]

    # Allocate data for the csc object
    csc_data = np.zeros(len(data), dtype=type(data[0]))

    # Scan through the CSR data structure
    for i in range(n):
        for jp in range(rowp[i], rowp[i + 1]):
            # Set the new row location in the CSC data structure
            j = cols[jp]
            csc_data[colp[j]] = data[jp]
            rows[colp[j]] = i
            colp[j] += 1

    # Reset the colp pointer
    for j in range(m, 0, -1):
        colp[j] = colp[j - 1]
    colp[0] = 0

    return {"csc": [colp, rows, csc_data], "shape": [n, m]}


def convertToDense(mat: Union[dict, spmatrix, ndarray]) -> ndarray:
    """
    Take a pyopsparse sparse matrix definition and convert back to a dense
    format. This is typically the final step for optimizers with dense constraint
    jacibians.

    Parameters
    ----------
    mat : dict
       A sparse matrix representation. Should be in CSR format for best
       efficiency

    Returns
    -------
    newMat : array
        A dense numpy array of the same matrix
    """

    mat = convertToCSR(mat)
    newMat = np.zeros(mat["shape"])
    data = mat["csr"][IDATA]
    colInd = mat["csr"][ICOLIND]
    rowp = mat["csr"][IROWP]
    for i in range(mat["shape"][0]):
        for j in range(rowp[i], rowp[i + 1]):
            newMat[i, colInd[j]] = data[j]
    return newMat


def scaleColumns(mat: dict, factor):
    """
    Scale the columns of the matrix. Must be CSR format
    """
    if not isinstance(mat, dict):
        raise TypeError("mat for scaleColumns must be pyoptsparse matrix format")
    if "csr" not in mat:
        raise ValueError("scaleColumns only works for CSR pyoptsparse matrix format")
    if mat["shape"][1] != len(factor):
        raise ValueError("Length of factor is incorrect")
    for i in range(mat["shape"][0]):
        iStart = mat["csr"][IROWP][i]
        iEnd = mat["csr"][IROWP][i + 1]
        mat["csr"][IDATA][iStart:iEnd] *= factor[mat["csr"][ICOLIND][iStart:iEnd]]


def scaleRows(mat: dict, factor):
    """
    Scale the rows of the matrix. Must be CSR format
    """
    if not isinstance(mat, dict):
        raise TypeError("mat for scaleRows must be pyoptsparse matrix format")
    if "csr" not in mat:
        raise ValueError("scaleRows only works for CSR pyoptsparse matrix format")
    if mat["shape"][0] != len(factor):
        raise ValueError("Length of factor is incorrect")
    for i in range(mat["shape"][0]):
        iStart = mat["csr"][IROWP][i]
        iEnd = mat["csr"][IROWP][i + 1]
        mat["csr"][IDATA][iStart:iEnd] *= factor[i]


def extractRows(mat: dict, indices):
    """
    Extract the rows defined by 'indices' and return
    a new CSR matrix.

    Parameters
    ----------
    mat : dict
        pyoptsparse matrix CSR format
    indices : list/array of integer
        The rows the user wants to extract

    Returns
    -------
    newMat : dic
       pyoptsparse CSR matrix
    """
    rowp = mat["csr"][IROWP]
    cols = mat["csr"][ICOLIND]
    data = mat["csr"][IDATA]
    m = mat["shape"][1]
    nn = len(indices)
    nrowp = np.zeros(nn + 1, "intc")

    # Count up the size of everything
    size = 0
    for i in range(nn):
        size += rowp[indices[i] + 1] - rowp[indices[i]]
        nrowp[i + 1] = size

    # Create the new columns and data arrays
    ncols = np.zeros(size, "intc")
    ndata = np.zeros(size, dtype=type(data[0]))

    # Re-indices the new columns
    for i in range(nn):
        ncols[nrowp[i] : nrowp[i + 1]] = cols[rowp[indices[i]] : rowp[indices[i] + 1]]
        ndata[nrowp[i] : nrowp[i + 1]] = data[rowp[indices[i]] : rowp[indices[i] + 1]]

    return {"csr": [nrowp, ncols, ndata], "shape": [nn, m]}


def _denseToCOO(arr: ndarray) -> dict:
    """
    Return a COO array that is a COO representation of the dense numpy
    array, arr

    Parameters
    ----------
    arr : numpy array

    Returns
    -------
    dict : mat
        The pyoptsparse representation of a sparse matrix
    """
    nRows = arr.shape[0]
    nCols = arr.shape[1]
    data = arr.flatten()
    cols = np.mod(np.arange(nRows * nCols), nCols)
    rows = np.arange(nRows * nCols) // nCols
    return {"coo": [rows, cols, data], "shape": [nRows, nCols]}


def _csr_to_coo(mat: dict) -> dict:
    """
    Convert the given CSR matrix to a COO format

    Parameters
    ----------
    mat : dict
       pyoptsparse matrix definition
    """

    # This is straight forward - just expand out the rows
    rowp = mat["csr"][IROWP]
    cols = mat["csr"][ICOLIND]
    data = mat["csr"][IDATA]
    coo_rows = np.zeros(len(cols), "intc")
    coo_cols = np.array(cols, "intc")

    for i in range(mat["shape"][0]):
        coo_rows[rowp[i] : rowp[i + 1]] = i

    coo_data = np.array(data)

    return {"coo": [coo_rows, coo_cols, coo_data], "shape": mat["shape"]}


def _csc_to_coo(mat: dict) -> dict:
    """
    Convert the given CSC matrix to a COO format

    Parameters
    ----------
    mat : dict
       pyoptsparse matrix definition
    """

    # This is straight forward - just expand out the rows
    colp = mat["csc"][ICOLP]
    rows = mat["csc"][IROWIND]
    data = mat["csc"][IDATA]

    # This is straight forward - just expand out the columns
    coo_rows = np.array(rows, "intc")
    coo_cols = np.zeros(len(rows), "intc")

    for j in range(mat["shape"][1]):
        coo_cols[colp[j] : colp[j + 1]] = j

    coo_data = np.array(data)

    return {"coo": [coo_rows, coo_cols, coo_data], "shape": mat["shape"]}


def _broadcast_to_array(name: str, value: ArrayType, n_values: int, allow_none: bool = False):
    """
    Broadcast an input to an array with a specified length

    Parameters
    ----------
    name : str
        The name of the input. This is only used in the error message emitted.
    value : float, list[float], numpy array
        The input value
    n_values : int
        The number of values
    allow_none : bool, optional
        Whether to allow `None` in the input/output, by default False

    Returns
    -------
    NDArray
        An array with the shape ``(n_values)``

    Raises
    ------
    ValueError
        If either the input is not broadcastable, or if the input contains None and ``allow_none=False``.

    Warnings
    --------
    Note that the default value for ``allow_none`` is False.
    """
    try:
        value = np.broadcast_to(value, n_values)
    except ValueError as e:
        raise ValueError(
            f"The '{name}' argument is invalid. It must be None, a scalar, or a list/array or length {n_values}."
        ) from e
    if not allow_none and any([i is None for i in value]):
        raise ValueError(f"The {name} argument cannot be 'None'.")
    return value


def try_import_compiled_module_from_path(
    module_name: str, path: Optional[str] = None, raise_warning: bool = False
) -> Union[types.ModuleType, str]:
    """
    Attempt to import a module from a given path.

    Parameters
    ----------
    module_name : str
        The name of the module
    path : Optional[str]
        The path to import from. If None, the default ``sys.path`` is used.
    raise_warning : bool
        If true, raise an import warning. By default false.

    Returns
    -------
    Union[types.ModuleType, str]
        If importable, the imported module is returned.
        If not importable, the error message is instead returned.
    """
    orig_path = sys.path
    if path is not None:
        path = os.path.abspath(os.path.expandvars(os.path.expanduser(path)))
        sys.path = [path]
    try:
        module = importlib.import_module(module_name)
    except ImportError as e:
        if raise_warning:
            warnings.warn(
                f"{module_name} module could not be imported from {path}.",
                stacklevel=2,
            )
        module = str(e)
    finally:
        sys.path = orig_path
    return module
