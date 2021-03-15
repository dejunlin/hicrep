#! /usr/bin/env python3
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright 2019 dejunlin <dejun.lin@gmail.com>
# Usage: hicrep.py
# Description: Compute HiCRep reproducibility stratum-corrected correlation score (SCCS).
# Reference: Genome Res. 2017 Nov;27(11):1939-1949. doi: 10.1101/gr.220640.117
# The algorithm first normalizes the input contact matrices by the total
# number of contacts and then for each chromosome: 1) mean-filter the input
# matrices with an input window size; 2) exclude common zero entries in
# the input matrices; 3) compute the SCC score. It doesn't have the
# procedure to bootstrap the window-size parameter
#
# Distributed under terms of the GNU General Public License v3.0.
import os
from deprecated import deprecated
import numpy as np
import scipy.sparse as sp
import math
import sys
import warnings
import cooler
from hicrep.utils import (
    readMcool, cool2pixels, getSubCoo,
    trimDiags, meanFilterSparse, varVstran,
    resample, upperDiagCsr, coolerInfo
    )

@deprecated("Use sccByDiag instead")
def sccOfDiag(diag1: np.ndarray, diag2: np.ndarray):
    """Get the correlation coefficient and weight of two input
    diagonal arrays

    Args:
        diag1: `np.ndarray` input array 1
        diag2: `np.ndarray` input array 2

    Returns:
       tuple of 2 floats, the Pearson's correlation rho and weight
    """
    # remove common zeros
    idxNZ = np.where((diag1 != 0.0) | (diag2 != 0.0))[0]
    iN = idxNZ.size
    if iN <= 2:
        return (np.nan, np.nan)
    iDiagNZ1 = diag1[idxNZ]
    iDiagNZ2 = diag2[idxNZ]
    rho = np.corrcoef(iDiagNZ1, iDiagNZ2)[0, 1]
    iDiagVarVstran1 = varVstran(iDiagNZ1.size)
    iDiagVarVstran2 = varVstran(iDiagNZ2.size)
    ws = iN * np.sqrt(iDiagVarVstran1 * iDiagVarVstran2)
    if math.isnan(rho) or math.isnan(ws):
        return (np.nan, np.nan)
    return (rho, ws)

def sccByDiag(m1: sp.coo_matrix, m2: sp.coo_matrix, nDiags: int):
    """Compute diagonal-wise hicrep SCC score for the two input matrices up to
    nDiags diagonals


    Args:
        m1 (sp.coo_matrix): input contact matrix 1
        m2 (sp.coo_matrix): input contact matrix 2
        nDiags (int): compute SCC scores for diagonals whose index is in the
        range of [1, nDiags)
    Returns: `float` hicrep SCC scores
    """
    # convert each diagonal to one row of a csr_matrix in order to compute
    # diagonal-wise correlation between m1 and m2
    m1D = upperDiagCsr(m1, nDiags)
    m2D = upperDiagCsr(m2, nDiags)
    nSamplesD = (m1D + m2D).getnnz(axis=1)
    rowSumM1D = m1D.sum(axis=1).A1
    rowSumM2D = m2D.sum(axis=1).A1
    # ignore zero-division warnings because the corresponding elements in the
    # output don't contribute to the SCC scores
    with np.errstate(divide='ignore', invalid='ignore'):
        cov = m1D.multiply(m2D).sum(axis=1).A1 - rowSumM1D * rowSumM2D / nSamplesD
        rhoD = cov / np.sqrt(
            (m1D.power(2).sum(axis=1).A1 - np.square(rowSumM1D) / nSamplesD ) *
            (m2D.power(2).sum(axis=1).A1 - np.square(rowSumM2D) / nSamplesD ))
        wsD = nSamplesD * varVstran(nSamplesD)
        # Convert NaN and Inf resulting from div by 0 to zeros.
        # posinf and neginf added to fix behavior seen in 4DN datasets
        # 4DNFIOQLTI9G and DNFIH7MQHOR at 5kb where inf would be reported
        # as an SCC score
        wsNan2Zero = np.nan_to_num(wsD, copy=True, posinf=0.0, neginf=0.0)
        rhoNan2Zero = np.nan_to_num(rhoD, copy=True, posinf=0.0, neginf=0.0)

    return rhoNan2Zero @ wsNan2Zero / wsNan2Zero.sum()


def hicrepSCC(cool1: cooler.api.Cooler, cool2: cooler.api.Cooler,
              h: int, dBPMax: int, bDownSample: bool):
    """Compute hicrep score between two input Cooler contact matrices

    Args:
        cool1: `cooler.api.Cooler` Input Cooler contact matrix 1
        cool2: `cooler.api.Cooler` Input Cooler contact matrix 2
        h: `int` Half-size of the mean filter used to smooth the
        input matrics
        dBPMax `int` Only include contacts that are at most this genomic
        distance (bp) away
        bDownSample: `bool` Down sample the input with more contacts
        to the same number of contacts as in the other input

    Returns:
        `float` scc scores for each chromosome
    """
    binSize1 = cool1.binsize
    binSize2 = cool2.binsize
    assert binSize1 == binSize2,\
        f"Input cool files have different bin sizes"
    assert coolerInfo(cool1, 'nbins') == coolerInfo(cool2, 'nbins'),\
        f"Input cool files have different number of bins"
    assert coolerInfo(cool1, 'nchroms') == coolerInfo(cool2, 'nchroms'),\
        f"Input cool files have different number of chromosomes"
    assert (cool1.chroms()[:] == cool2.chroms()[:]).all()[0],\
        f"Input file have different chromosome names"
    binSize = binSize1
    bins1 = cool1.bins()
    bins2 = cool2.bins()
    if binSize is None:
        # sometimes bin size can be None, e.g., input cool file has
        # non-uniform size bins.
        assert np.all(bins1[:] == bins2[:]),\
            f"Input cooler files don't have a unique bin size most likely "\
            f"because non-uniform bin size was used and the bins are defined "\
            f"differently for the two input cooler files"
        # In that case, use the median bin size
        binSize = int(np.median((bins1[:]["end"] - bins1[:]["start"]).values))
        warnings.warn(f"Input cooler files don't have a unique bin size most "\
                      f"likely because non-uniform bin size was used. HicRep "\
                      f"will use median bin size from the first cooler file "\
                      f"to determine maximal diagonal index to include", RuntimeWarning)
    if dBPMax == -1:
        # this is the exclusive upper bound
        dMax = coolerInfo(cool1, 'nbins')
    else:
        dMax = dBPMax // binSize + 1
    assert dMax > 1, f"Input dBPmax is smaller than binSize"
    p1 = cool2pixels(cool1)
    p2 = cool2pixels(cool2)
    # get the total number of contacts as normalizing constant
    n1 = coolerInfo(cool1, 'sum')
    n2 = coolerInfo(cool2, 'sum')
    chrNames = cool1.chroms()[:]['name'].to_numpy()
    # filter out mitochondria chromosome
    chrNames = np.array([name for name in chrNames if name != 'M'])
    scc = np.full(chrNames.shape[0], -2.0)
    for iChr in range(chrNames.shape[0]):
        chrName = chrNames[iChr]
        # normalize by total number of contacts
        mS1 = getSubCoo(p1, bins1, chrName)
        assert mS1.size > 0, "Contact matrix 1 of chromosome %s is empty" % (chrName)
        assert mS1.shape[0] == mS1.shape[1],\
            "Contact matrix 1 of chromosome %s is not square" % (chrName)
        mS2 = getSubCoo(p2, bins2, chrName)
        assert mS2.size > 0, "Contact matrix 2 of chromosome %s is empty" % (chrName)
        assert mS2.shape[0] == mS2.shape[1],\
            "Contact matrix 2 of chromosome %s is not square" % (chrName)
        assert mS1.shape == mS2.shape,\
            "Contact matrices of chromosome %s have different input shape" % (chrName)
        nDiags = mS1.shape[0] if dMax < 0 else min(dMax, mS1.shape[0])
        rho = np.full(nDiags, np.nan)
        ws = np.full(nDiags, np.nan)
        # remove major diagonal and all the diagonals >= nDiags
        # to save computation time
        m1 = trimDiags(mS1, nDiags, False)
        m2 = trimDiags(mS2, nDiags, False)
        del mS1
        del mS2
        if bDownSample:
            # do downsampling
            size1 = m1.sum()
            size2 = m2.sum()
            if size1 > size2:
                m1 = resample(m1, size2).astype(float)
            elif size2 > size1:
                m2 = resample(m2, size1).astype(float)
        else:
            # just normalize by total contacts
            m1 = m1.astype(float) / n1
            m2 = m2.astype(float) / n2
        if h > 0:
            # apply smoothing
            m1 = meanFilterSparse(m1, h)
            m2 = meanFilterSparse(m2, h)
        scc[iChr] = sccByDiag(m1, m2, nDiags)
    return scc
