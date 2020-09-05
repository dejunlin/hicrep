#! /usr/bin/env python3
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright 2019 dejunlin <dejun.lin@gmail.com>
# Usage: utils.py
# Description: Utility functions
#
# Distributed under terms of the MIT license.
import numpy as np
import pandas as pd
import cooler
import h5py
import math
import scipy.sparse as sp

def readMcool(fmcool: str, binSize: int):
    """Read from a mcool or cool file and return the Cooler object

    Args:
        fmcool: Input file name
        binSize: Bin size to select from the mcool file. If this value
        is <= 0, the input will be treated as a cool file instead

    Returns:
        cooler.api.Cooler object
    """
    mcool = h5py.File(fmcool, 'r')
    if binSize > 0:
        return cooler.Cooler(mcool['resolutions'][str(binSize)]), binSize
    else:
        cool = cooler.Cooler(mcool)
        return cool, cool.binsize


def cool2pixels(cool: cooler.api.Cooler):
    """Return the contact matrix in "pixels" format

    Args:
        cool: Input cooler object

    Returns:
        cooler.core.RangeSelector2D object
    """
    return cool.matrix(as_pixels=True, balance=False, sparse=True)


def pixels2Coo(df: pd.DataFrame, bins: pd.DataFrame):
    """Convert Cooler's contact matrix in "pixels" DataFrame to
    scipy coo_matrix. The "pixels" format is a 3-column DataFrame:
    'bin1_id', 'bin2_id', 'counts' for each unique contact

    Args:
        df: Input DataFrame
        bins: Cooler bins for the contacts

    Returns:
        coo_matrix of the input
    """
    binOffset = bins.index[0]
    nBins = bins.shape[0]
    df['bin1_id'] -= binOffset
    df['bin2_id'] -= binOffset
    return sp.coo_matrix((df['count'].to_numpy(),
                          (df['bin1_id'].to_numpy(), df['bin2_id'].to_numpy())),
                         shape=(nBins, nBins))


def getSubCoo(pixels: cooler.core.RangeSelector2D, bins: cooler.core.RangeSelector1D,
              regionStr: str):
    """Fetch a region from Cooler a contact matrix and return it as a
    coo_matrix

    Args:
        pixels: Input Cooler range selector object of the contact matrix
        pixels: Input Cooler range selector object of the bin definition
        regionStr: String for selecting genomic region

    Returns:
        coo_matrix contact matrix corresponding to the input region
    """
    mSub = pixels.fetch(regionStr)
    # Assume Cooler always use upper triangle
    assert (mSub['bin1_id'] <= mSub['bin2_id']).all(),\
        f"Contact matrix of region {regionStr} has lower-triangle entries"
    binsSub = bins.fetch(regionStr)
    return pixels2Coo(mSub, binsSub)


def trimDiags(a: sp.coo_matrix, iDiagMax: int, bKeepMain: bool):
    """Remove diagonal elements whose diagonal index is >= iDiagMax
    or is == 0

    Args:
        a: Input scipy coo_matrix
        iDiagMax: Diagonal offset cutoff
        bKeepMain: If true, keep the elements in the main diagonal;
        otherwise remove them

    Returns:
        coo_matrix with the specified diagonals removed
    """
    gDist = np.abs(a.row - a.col)
    idx = np.where((gDist < iDiagMax) & (bKeepMain | (gDist != 0)))
    return sp.coo_matrix((a.data[idx], (a.row[idx], a.col[idx])),
                         shape=a.shape, dtype=a.dtype)


def meanFilterSparse(a: sp.coo_matrix, h: int):
    """Apply a mean filter to an input sparse matrix. This convolves
    the input with a kernel of size 2*h + 1 with constant entries and
    subsequently reshape the output to be of the same shape as input

    Args:
        a: `sp.coo_matrix`, Input matrix to be filtered
        h: `int` half-size of the filter

    Returns:
        `sp.coo_matrix` filterd matrix
    """
    assert h > 0, "meanFilterSparse half-size must be greater than 0"
    assert sp.issparse(a) and a.getformat() == 'coo',\
        "meanFilterSparse input matrix is not scipy.sparse.coo_matrix"
    assert a.shape[0] == a.shape[1],\
        "meanFilterSparse cannot handle non-square matrix"
    fSize = 2 * h + 1
    # filter is a square matrix of constant 1 of shape (fSize, fSize)
    shapeOut = np.array(a.shape) + fSize - 1
    mToeplitz = sp.diags(np.ones(fSize),
                         np.arange(-fSize+1, 1),
                         shape=(shapeOut[1], a.shape[1]),
                         format='csr')
    ans = sp.coo_matrix((mToeplitz @ a) @ mToeplitz.T)
    # remove the edges since we don't care about them if we are smoothing
    # the matrix itself
    ansNoEdge = ans.tocsr()[h:(h+a.shape[0]), h:(h+a.shape[1])].tocoo()
    # Assign different number of neighbors to the edge to better
    # match what the original R implementation of HiCRep does
    rowDist2Edge = np.minimum(ansNoEdge.row, ansNoEdge.shape[0] - 1 - ansNoEdge.row)
    nDim1 = h + 1 + np.minimum(rowDist2Edge, h)
    colDist2Edge = np.minimum(ansNoEdge.col, ansNoEdge.shape[1] - 1 - ansNoEdge.col)
    nDim2 = h + 1 + np.minimum(colDist2Edge, h)
    nNeighbors = nDim1 * nDim2
    ansNoEdge.data /= nNeighbors
    return ansNoEdge

def varVstran(a: np.ndarray):
    """
    Calculate the variance of variance-stabilizing transformed
    (or `vstran()` in the original R implementation) data. The `vstran()` turns
    the input data into ranks, whose variance is only a function of the input
    size:
        ```
        var(1/n, 2/n, ..., n/n) = (1 - 1/(n^2))/12
        ```
    or with Bessel's correction:
        ```
        var(1/n, 2/n, ..., n/n, ddof=1) = (1 + 1.0/n)/12
        ```
    See section "Variance stabilized weights" in reference for more detail:
    https://genome.cshlp.org/content/early/2017/10/06/gr.220640.117

    Args:
        a (np.ndarray): Input data
    Returns: `float` variance of the ranked input data with Bessel's correction
    """
    return (1 + 1.0 / a.shape[0]) / 12.0


def resample(m: sp.coo_matrix, size: int):
    """Resample with replacement the input matrix so that the
    resulting matrix sum to the given size
    Args:
        m: `sp.coo_matrix` Input matrix
        size: Resulting matrix sum to this number

    Returns:
        resampled matrix
    """
    bins = np.arange(m.data.size)
    p = m.data / m.data.sum()
    samples = np.random.choice(bins, size=size, p=p)
    sampledData = np.bincount(samples, minlength=bins.size)
    ans = sp.coo_matrix((sampledData, (m.row, m.col)), shape=m.shape)
    ans.eliminate_zeros()
    return ans
