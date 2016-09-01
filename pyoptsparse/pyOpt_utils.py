#!/usr/bin/env python
'''
pyOptSparse_utils

Holds a minimal set of sparse-matrix type routines for pyOptSparse. This
is designed to replace the HORRENDOUS scipy sparse matrix format. The
with scipy.sparse is that is the NO way to enforce a constant sparsity
structure which is required for the optimizers. We use a very simple
dictionary format to represent the three most common forms of sparse matrices:

mat = {'coo':[row, col, data], 'shape':[nrow, ncols]} # A coo matrix
mat = {'csr':[rowp, colind, data], 'shape':[nrow, ncols]} # A csr matrix
mat = {'coo':[colp, rowind, data], 'shape':[nrow, ncols]} # A csc matrix

Copyright (c) 2008-2013
All rights reserved.

Developers:
-----------
- Dr. Gaetan K.W. Kenway (GKK)
- Dr. Graeme J. Kennedy (GJK)

History
-------
    v. 1.0  - Initial Class Creation (GKK, 2014)
'''
import numpy
import warnings
from .pyOpt_error import Error
# Define index memonics
IROW = 0
ICOL = 1

IROWP = 0
ICOLIND = 1

ICOLP = 0
IROWIND = 1

IDATA = 2
__all__ = ['convertToCOO', 'convertToCSR', 'convertToCSC', 'convertToDense',
           'multCOO',
           'scaleColumns', 'scaleRows', 'extractRows', 'IROW', 'ICOL',
           'IROWP', 'ICOLIND', 'ICOLP', 'IROWIND', 'IDATA']

def convertToCOO(mat):
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
        if 'coo' in mat:
            return mat
        if 'csr' in mat:
            return _csr_to_coo(mat)
        elif 'csc' in mat:
            return _csc_to_coo(mat)
    else:
        # Try to do it with a scipy sparse matrix:
        try:
            from scipy import sparse
            if sparse.issparse(mat):
                warnings.warn("Using scipy.sparse matrices with pyOptSparse "
                              "in VERY STRONGLY discouraged. Please use the "
                              "simplified pyoptsparse format which allows for "
                              "fixed sparsity structure and explict zeros in "
                              "the matrix. There is no way to guarantee "
                              "a fixed sparsity structure with scipy "
                              "matrices which is what the underlying "
                              "optimizers require. Using scipy.sparse "
                              "matrices may cause unexpected errors.")

                mat = mat.tocoo()
                return {'coo':[mat.row, mat.col, mat.data], 'shape':mat.shape}
        except:
            pass

        # Now try to do it with a numpy matrix:
        try:
            return _denseToCOO(numpy.atleast_2d(numpy.array(mat)))
        except:
            raise Error("Unknown matrix format. Must be a dense numpy "
                        "array or a pyoptsparse sparce matrix format of "
                        "COO, CSR or CSC. See documentation for correct "
                        "format. Supplied Matrix is: %s"% repr(mat))

def convertToCSR(mat):
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

    if 'csr' in mat:
        return mat

    mat = convertToCOO(mat)
    n = mat['shape'][0]
    m = mat['shape'][1]
    rows = mat['coo'][IROW]
    cols = mat['coo'][ICOL]
    data = mat['coo'][IDATA]

    rowp = numpy.zeros(n+1, dtype='intc')

    # Count up the number of times things are index
    for row in rows:
        rowp[row+1] += 1

    # Set up the array as a pointer
    for i in range(1, n+1):
        rowp[i] += rowp[i-1]

    ncols = numpy.zeros(rowp[-1], dtype='intc')
    ndata = numpy.zeros(rowp[-1], dtype=type(data[0]))

    # Now, add all the values and the data
    for i in range(len(rows)):
        r = rows[i]
        ncols[rowp[r]] = cols[i]
        ndata[rowp[r]] = data[i]
        rowp[r] += 1

    # Readjust the pointer
    for i in range(n, 0, -1):
        rowp[i] = rowp[i-1]
    rowp[0] = 0

    return {'csr':[rowp, ncols, ndata], 'shape':[n, m]}

def convertToCSC(mat):
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
    if 'csc' in mat:
        return mat

    mat = convertToCSR(mat)
    n = mat['shape'][0]
    m = mat['shape'][1]
    rowp = mat['csr'][IROWP]
    cols = mat['csr'][ICOLIND]
    data = mat['csr'][IDATA]

    # Allocate the new arrays
    colp = numpy.zeros(m+1, 'intc')
    rows = numpy.zeros(len(cols), 'intc')

    # Count up the number of references to each column
    for col in cols:
        colp[col+1] += 1

    # Set colp so that it is now a pointer
    for i in range(1, m):
        colp[i] += colp[i-1]

    # Allocate data for the csc object
    csc_data = numpy.zeros(len(data), dtype=type(data[0]))

    # Scan through the CSR data structure
    for i in range(n):
        for jp in range(rowp[i], rowp[i+1]):
            # Set the new row location in the CSC data structure
            j = cols[jp]
            csc_data[colp[j]] = data[jp]
            rows[colp[j]] = i
            colp[j] += 1

    # Reset the colp pointer
    for j in range(m, 0, -1):
        colp[j] = colp[j-1]
    colp[0] = 0

    return {'csc':[colp, rows, csc_data], 'shape':[n, m]}

def convertToDense(mat):
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
    newMat = numpy.zeros((mat['shape']))
    data = mat['csr'][IDATA]
    colInd = mat['csr'][ICOLIND]
    rowp = mat['csr'][IROWP]
    for i in range(mat['shape'][0]):
        for j in range(rowp[i], rowp[i+1]):
            newMat[i, colInd[j]] = data[j]
    return newMat

def scaleColumns(mat, factor):
    """ d=
    Scale the columns of the matrix. Must be CSR format
    """
    if not isinstance(mat, dict):
        raise Error("mat for scaleColumbs must be pyoptsparse matrix format")
    if 'csr' not in mat:
        raise Error("scaleColumns only works for CSR pyoptsparse matrix format")
    if mat['shape'][1] != len(factor):
        raise Error("Length of factor is incorrect")
    for i in range(mat['shape'][0]):
        iStart = mat['csr'][IROWP][i]
        iEnd = mat['csr'][IROWP][i+1]
        for j in range(iStart, iEnd):
            mat['csr'][IDATA][j] *= factor[mat['csr'][ICOLIND][j]]

def scaleRows(mat, factor):
    """
    Scale the rows of the matrix. Must be CSR format
    """
    if not isinstance(mat, dict):
        raise Error("mat for scaleRows must be pyoptsparse matrix format")
    if 'csr' not in mat:
        raise Error("scaleRows only works for CSR pyoptsparse matrix format")
    if mat['shape'][0] != len(factor):
        raise Error("Length of factor is incorrect")
    for i in range(mat['shape'][0]):
        iStart = mat['csr'][IROWP][i]
        iEnd = mat['csr'][IROWP][i+1]
        mat['csr'][IDATA][iStart:iEnd] *= factor[i]

def extractRows(mat, indices):
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
    mat = convertToCSR(mat)
    rowp = mat['csr'][IROWP]
    cols = mat['csr'][ICOLIND]
    data = mat['csr'][IDATA]
    m = mat['shape'][1]
    nn = len(indices)
    nrowp = numpy.zeros(nn+1, 'intc')

    # Count up the size of everything
    size = 0
    for i in range(nn):
        size += rowp[indices[i]+1] - rowp[indices[i]]
        nrowp[i+1] = size

    # Create the new columns and data arrays
    ncols = numpy.zeros(size, 'intc')
    ndata = numpy.zeros(size, dtype=type(data[0]))

    # Re-indices the new columns
    for i in range(nn):
        ncols[nrowp[i]:nrowp[i+1]] = cols[rowp[indices[i]]:rowp[indices[i]+1]]
        ndata[nrowp[i]:nrowp[i+1]] = data[rowp[indices[i]]:rowp[indices[i]+1]]

    return {'csr':[nrowp, ncols, ndata], 'shape':[nn, m]}

def multCOO(mat, x):
    """
    Inefficient matrix-vector multiply for COO matrix

    Parameters
    ----------
    mat : dict
        pyoptsparse coo mat
    x : numpy array
        vector to multiply with

    Returns
    -------
    y : numpy array
        resulting vector
        """
    y = numpy.zeros(mat['shape'][0], dtype=type(x[0]))
    rows = mat['coo'][IROW]
    cols = mat['coo'][ICOL]
    data = mat['coo'][IDATA]
    for i in range(len(data)):
        y[rows[i]] += data[i]*x[cols[i]]

    return y

def _denseToCOO(arr):
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
    cols = numpy.mod(numpy.arange(nRows*nCols), nCols)
    rows = numpy.arange(nRows*nCols)//nCols
    return {'coo':[rows, cols, data], 'shape':[nRows, nCols]}

def _csr_to_coo(mat):
    """
    Convert the given CSR matrix to a COO format

    Parameters
    ----------
    mat : dict
       pyoptsparse matrix definition
    """

    # This is straight forward - just expand out the rows
    rowp = mat['csr'][IROWP]
    cols = mat['csr'][ICOLIND]
    data = mat['csr'][IDATA]
    coo_rows = numpy.zeros(len(cols), 'intc')
    coo_cols = numpy.array(cols, 'intc')

    for i in range(mat['shape'][0]):
        coo_rows[rowp[i]:rowp[i+1]] = i

    coo_data = numpy.array(data)

    return {'coo':[coo_rows, coo_cols, coo_data], 'shape':mat['shape']}

def _csc_to_coo(mat):
    """
    Convert the given CSC matrix to a COO format

    Parameters
    ----------
    mat : dict
       pyoptsparse matrix definition
    """

    # This is straight forward - just expand out the rows
    colp = mat['csc'][ICOLP]
    rows = mat['csc'][IROWIND]
    data = mat['csc'][IDATA]

    # This is straight forward - just expand out the columns
    coo_rows = numpy.array(rows, 'intc')
    coo_cols = numpy.zeros(len(rows), 'intc')

    for j in range(mat['shape'][1]):
        coo_cols[colp[j]:colp[j+1]] = j

    coo_data = numpy.array(data)

    return {'coo':[coo_rows, coo_cols, coo_data], 'shape':mat['shape']}
