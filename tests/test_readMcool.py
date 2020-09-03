#! /usr/bin/env python3
#
# Testing for the helper function readMcool in hicrep/utils
#
import pytest
from hicrep.utils import readMcool

def testReadMcool():
    fmcool = "tests/data/human_hi-c/4DNFICQK4N8B.mcool"

    # Test that readMcool throws an exception if you try to read in an .mcool
    # file at an unavailable resolution
    with pytest.raises(KeyError):
        readMcool(fmcool, 12345)

    # Test that the read in .mcool file has the specified resolution and the
    # correct chromosome labels
    cool, binSize = readMcool(fmcool, 25000)
    assert cool.binsize == 25000,\
        f"Cooler object .mcool file {fmcool} has different resolution than requested"

    fmcool = "tests/data/human_hi-c/4DNFIQ5XCHDB.cool"

    # Test that readMcool throws an error if you try and read in a .cool file at
    # a different resolution than it's builtin
    with pytest.raises(KeyError):
        readMcool(fmcool, 25000)

    # Test that the read in .cool file has the correct built in binsize and the
    # correct chromosome labels
    cool, binSize = readMcool(fmcool, -1)
    assert cool.binsize == 500000,\
        f"Cooler object .mcool file {fmcool} has different resolution than requested"
