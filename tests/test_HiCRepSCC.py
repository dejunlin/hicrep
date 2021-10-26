#! /usr/bin/env python3
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright 2020 dejunlin <dejun.lin@gmail.com>
# Created at 2020-02-12 17:37 on dejunlin@threonine.gs.washington.edu
# Usage: testHiCRepSCC.py
# Description: Test HiCRep using public data
#
# Distributed under terms of the GNU General Public License v3.0.
import numpy as np
from hicrep.utils import (
    readMcool, cool2pixels, getSubCoo,
    trimDiags, meanFilterSparse,
    resample, coolerInfo
    )
from hicrep.hicrep import (
    sccOfDiag, hicrepSCC
    )

def testHumanHiC():
    fmcool1 = "tests/data/human_hi-c/4DNFITKCX2DO.cool"
    fmcool2 = "tests/data/human_hi-c/4DNFIQ5XCHDB.cool"
    binSize = 500000
    h = 0
    dBPMax = 5000000
    bDownSample = False
    cool1, _ = readMcool(fmcool1, -1)
    cool2, _ = readMcool(fmcool2, -1)

    #Test that the scc scores between a matrix and itself are always 1
    results = hicrepSCC(cool1, cool1, h, dBPMax, bDownSample)
    assert np.isclose(results, 1).all(),\
        f"SCC scores between {fmcool1} and itself are not 1"

    #Test that the scc scores agree with the previous R implementation
    results = hicrepSCC(cool1, cool2, h, dBPMax, bDownSample)
    #Values given by R implementation of hicrep with same parameters
    expected = np.array([0.73050741, 0.67516601, 0.65743544, 0.74419469,
                        0.75864553, 0.75172288, 0.79556107, 0.66194007,
                        0.7141874,  0.78722237, 0.77622226, 0.77451858,
                        0.73061994, 0.70921468, 0.7447885,  0.75176337,
                        0.77104526, 0.83237602, 0.79166534, 0.80379132,
                        0.7504225,  0.64200014, 0.84293773, 0.79261671])
    assert np.isclose(results, expected).all(),\
        f"SCC scores between {fmcool1} and {fmcool2} differ from those given by the R implementation"

    # Test the computation of a subset of chromosomes give the same results as
    # the whole set
    chrNames = ['chr1', 'chr2', 'chrX']
    resultsSub = hicrepSCC(cool1, cool2, h, dBPMax, bDownSample, chrNames)
    chrNamesAll = cool1.chroms()[:]['name'].tolist()
    chrNamesAll = [name for name in chrNamesAll if name != 'M']
    iChrs = np.where(np.isin(chrNamesAll, chrNames))[0]
    assert (results[iChrs] == resultsSub).all(), f"""
        SCC scores between {fmcool1} and {fmcool2} on chromosome subset
        {chrNames} differ from those computed from the whole set. The whole
        genome results are: {results} and the subset indices are {iChrs}.
        """

    # Test the computation when excluding chromosomes give the same results as
    # the whole set
    exclNames = set(['chr1', 'chr2', 'chrX'])
    resultsExcl = hicrepSCC(cool1, cool2, h, dBPMax, bDownSample,
                            chrNamesAll, exclNames)
    chrNamesRemain = [name for name in chrNamesAll if name not in exclNames]
    ieChrs = np.where(np.isin(chrNamesAll, chrNamesRemain))[0]
    assert (results[ieChrs] == resultsExcl).all(), f"""
        SCC scores between {fmcool1} and {fmcool2} differ from those computed
        from the whole set when {chrNamesRemain} are excluded. The whole genome
        results are: {results} and the subset indices are {iChrs}.
        """


def testFlyHiC():
    fmcool1 = "tests/data/fly_hi-c/4DNFI8DRD739_bin100kb.cool"
    fmcool2 = "tests/data/fly_hi-c/4DNFIZ1ZVXC8_bin100kb.cool"
    binSize = -1
    h = 1
    dBPMax = 500000
    bDownSample = False
    cool1, binSize1 = readMcool(fmcool1, binSize)
    cool2, binSize2 = readMcool(fmcool2, binSize)
    assert coolerInfo(cool1, 'nbins') == coolerInfo(cool2, 'nbins'),\
        f"Input cool files {fmcool1} and {fmcool2} have different number of bins"
    assert binSize1 == binSize2,\
        f"Input cool files {fmcool1} and {fmcool2} have different bin sizes"
    assert coolerInfo(cool1, 'nchroms') == coolerInfo(cool2, 'nchroms'),\
        f"Input cool files {fmcool1} and {fmcool2} have different number of chromosomes"
    assert (cool1.chroms()[:] == cool2.chroms()[:]).all()[0],\
        f"Input file {fmcool1} and {fmcool2} have different chromosome names"
    results = hicrepSCC(cool1, cool2, h, dBPMax, bDownSample)
    expected = np.array([9.936753824600870e-01, 9.950138992224218e-01,
                         9.951519844417879e-01, 9.935973973292749e-01,
                         9.933660605077106e-01, 9.927681695925705e-01,
                         6.238132870270471e-01])
    assert np.isclose(results, expected).all()

    # Test the computation of a subset of chromosomes give the same results as
    # the whole set
    chrNames = ['chr2L', 'chr2R', 'chrX']
    resultsSub = hicrepSCC(cool1, cool2, h, dBPMax, bDownSample, chrNames)
    chrNamesAll = cool1.chroms()[:]['name'].tolist()
    chrNamesAll = [name for name in chrNamesAll if name != 'M']
    iChrs = np.where(np.isin(chrNamesAll, chrNames))[0]
    assert (results[iChrs] == resultsSub).all(), f"""
        SCC scores between {fmcool1} and {fmcool2} on chromosome subset
        {chrNames} differ from those computed from the whole set. The whole
        genome results are: {results} and the subset indices are {iChrs}.
        """

    # Test the computation when excluding chromosomes give the same results as
    # the whole set
    exclNames = set(['chr2L', 'chr2R', 'chrX'])
    resultsExcl = hicrepSCC(cool1, cool2, h, dBPMax, bDownSample,
                            chrNamesAll, exclNames)
    chrNamesRemain = [name for name in chrNamesAll if name not in exclNames]
    ieChrs = np.where(np.isin(chrNamesAll, chrNamesRemain))[0]
    assert (results[ieChrs] == resultsExcl).all(), f"""
        SCC scores between {fmcool1} and {fmcool2} differ from those computed
        from the whole set when {chrNamesRemain} are excluded. The whole genome
        results are: {results} and the subset indices are {iChrs}.
        """
