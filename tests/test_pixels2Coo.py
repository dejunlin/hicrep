#! /usr/bin/env python3
#
# Testing for the helper function readMcool in hicrep/utils
#
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
    sparse = pixels2Coo(pix, bins)

    #Check that the shape of the sparse matrix is n * n where n is the number of bins in the original cooler (at the specified resolution)
    assert sparse.shape == (cool.info['nbins'], cool.info['nbins']),\
        f"pixels2coo() returns a matrix of different dhape than the original cooler"

    #Check that the number of data points in the sparse matrix gotten from pixels2coo is equal to
    #the number specified in the original cooler file
    assert sparse.nnz == cool.info['nnz'],\
        f"pixels2coo() returns a matrix with a different number of nonzero entries then specified in the original cooler"
