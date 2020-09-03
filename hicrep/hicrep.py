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
# Distributed under terms of the MIT license.
import os
import numpy as np
import math
import sys
import warnings
import cooler
from hicrep.utils import (
    readMcool, cool2pixels, getSubCoo,
    trimDiags, meanFilterSparse, varVstran,
    resample
    )


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
    iDiagVarVstran1 = varVstran(iDiagNZ1)
    iDiagVarVstran2 = varVstran(iDiagNZ2)
    ws = iN * np.sqrt(iDiagVarVstran1 * iDiagVarVstran2)
    if math.isnan(rho) or math.isnan(ws):
        return (np.nan, np.nan)
    return (rho, ws)


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
    binSize1 = cool1.info['bin-size']
    binSize2 = cool2.info['bin-size']
    assert binSize1 == binSize2,\
        f"Input cool files have different bin sizes"
    assert cool1.info['nbins'] == cool2.info['nbins'],\
        f"Input cool files have different number of bins"
    assert cool1.info['nchroms'] == cool2.info['nchroms'],\
        f"Input cool files have different number of chromosomes"
    assert (cool1.chroms()[:] == cool2.chroms()[:]).all()[0],\
        f"Input file have different chromosome names"
    binSize = binSize1
    if dBPMax == -1:
        # In general, don't use the entire contact matrix because usually the last
        # few diagonals have very few valid data in it for computing Pearson's correlation
        warnings.warn(f"Using dBPMax == -1 risk numerical instability at farthest "\
                      f"diagonals for computing Pearon's correlation", RuntimeWarning)
        # this is the exclusive upper bound
        dMax = cool1.info['nbins']
    else:
        dMax = dBPMax // binSize + 1
    assert dMax > 1, f"Input dBPmax is smaller than binSize"
    p1 = cool2pixels(cool1)
    p2 = cool2pixels(cool2)
    bins1 = cool1.bins()
    bins2 = cool2.bins()
    # get the total number of contacts as normalizing constant
    n1 = p1[:]['count'].sum()
    n2 = p2[:]['count'].sum()
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
        # ignore the main diagonal iD == 0
        for iD in range(1, nDiags):
            iDiag1 = m1.diagonal(iD)
            iDiag2 = m2.diagonal(iD)
            rho[iD], ws[iD] = sccOfDiag(iDiag1, iDiag2)
        wsNan2Zero = np.nan_to_num(ws, copy=True)
        rhoNan2Zero = np.nan_to_num(rho, copy=True)
        scc[iChr] = rhoNan2Zero @ wsNan2Zero / wsNan2Zero.sum()
    return scc
