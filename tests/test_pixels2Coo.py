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
import pytest
from hicrep.utils import (
    readMcool, cool2pixels, pixels2Coo, getSubCoo
    )

def testPixels2Coo():
    fmcool = "tests/data/human_hi-c/4DNFICQK4N8B.mcool"
    cool, binSize = readMcool(fmcool, 25000)

    #cool2pixels gets implicitly tested here
    pix = cool2pixels(cool)[:]
    bins = cool.bins()[:]
    mSparse = pixels2Coo(pix, bins)

    # Check that the shape of the sparse matrix is n * n where n is the number
    # of bins in the original cooler (at the specified resolution)
    assert mSparse.shape == (cool.info['nbins'], cool.info['nbins']),\
        f"pixels2coo() returns a matrix of different dhape than the original cooler"

    # Check that the number of data points in the sparse matrix gotten from
    # pixels2coo is equal to #the number specified in the original cooler file
    assert mSparse.nnz == cool.info['nnz'],\
        f"pixels2coo() returns a matrix with a different number of "\
        f"nonzero entries then specified in the original cooler"
