#! /usr/bin/env python3
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright 2021 dejun <dejun.lin@gmail.com>
# Created at 2021-03-15 12:28 on dejun@dejun-GS60-2QE
# Usage: test_coolerInfo.py
# Description: Test the util function coolerInfo()
#
# Distributed under terms of the GNU General Public License v3.0.
from hicrep.utils import (
    readMcool, coolerInfo
    )


def testHumanHiCInfo():
    fmcool1 = "tests/data/human_hi-c/4DNFITKCX2DO.cool"
    fmcool2 = "tests/data/human_hi-c/4DNFIQ5XCHDB.cool"
    cool1, binSize1 = readMcool(fmcool1, -1)
    cool2, binSize2 = readMcool(fmcool2, -1)
    
    #Check various .info() fields for consistency
    assert coolerInfo(cool1, 'bin-size') == binSize1,\
        f'coolerInfo() failed to retrieve metadata \'bin-size\' from {fmcool1}'
    assert coolerInfo(cool2, 'bin-size') == binSize2,\
        f'coolerInfo() failed to retrieve metadata \'bin-size\' from {fmcool2}'

    assert coolerInfo(cool1, 'sum') == cool1.pixels()['count'][:].sum(),\
        f'coolerInfo() failed to retrieve metadata \'sum\' from {fmcool1}'
    assert coolerInfo(cool2, 'sum') == cool2.pixels()['count'][:].sum(),\
        f'coolerInfo() failed to retrieve metadata \'sum\' from {fmcool2}'

    assert coolerInfo(cool1, 'nbins') == cool1.bins().shape[0],\
        f'coolerInfo() failed to retrieve metadata \'nbins\' from {fmcool1}'
    assert coolerInfo(cool2, 'nbins') == cool2.bins().shape[0],\
        f'coolerInfo() failed to retrieve metadata \'nbins\' from {fmcool2}'

    assert coolerInfo(cool1, 'nnz') == cool1.pixels().shape[0],\
        f'coolerInfo() failed to retrieve metadata \'nnz\' from {fmcool1}'
    assert coolerInfo(cool2, 'nnz') == cool2.pixels().shape[0],\
        f'coolerInfo() failed to retrieve metadata \'nnz\' from {fmcool2}'

    assert coolerInfo(cool1, 'nchroms') == cool1.chroms().shape[0],\
        f'coolerInfo() failed to retrieve metadata \'nchroms\' from {fmcool1}'
    assert coolerInfo(cool2, 'nchroms') == cool2.chroms().shape[0],\
        f'coolerInfo() failed to retrieve metadata \'nchroms\' from {fmcool2}'


def testFlyHiC():
    fmcool1 = "tests/data/fly_hi-c/4DNFI8DRD739_bin100kb.cool"
    fmcool2 = "tests/data/fly_hi-c/4DNFIZ1ZVXC8_bin100kb.cool"
    cool1, binSize1 = readMcool(fmcool1, -1)
    cool2, binSize2 = readMcool(fmcool2, -1)

    #Check various .info() fields for consistency
    assert coolerInfo(cool1, 'bin-size') == binSize1,\
        f'coolerInfo() failed to retrieve metadata \'bin-size\' from {fmcool1}'
    assert coolerInfo(cool2, 'bin-size') == binSize2,\
        f'coolerInfo() failed to retrieve metadata \'bin-size\' from {fmcool2}'

    assert coolerInfo(cool1, 'sum') == cool1.pixels()['count'][:].sum(),\
        f'coolerInfo() failed to retrieve metadata \'sum\' from {fmcool1}'
    assert coolerInfo(cool2, 'sum') == cool2.pixels()['count'][:].sum(),\
        f'coolerInfo() failed to retrieve metadata \'sum\' from {fmcool2}'

    assert coolerInfo(cool1, 'nbins') == cool1.bins().shape[0],\
        f'coolerInfo() failed to retrieve metadata \'nbins\' from {fmcool1}'
    assert coolerInfo(cool2, 'nbins') == cool2.bins().shape[0],\
        f'coolerInfo() failed to retrieve metadata \'nbins\' from {fmcool2}'

    assert coolerInfo(cool1, 'nnz') == cool1.pixels().shape[0],\
        f'coolerInfo() failed to retrieve metadata \'nnz\' from {fmcool1}'
    assert coolerInfo(cool2, 'nnz') == cool2.pixels().shape[0],\
        f'coolerInfo() failed to retrieve metadata \'nnz\' from {fmcool2}'

    assert coolerInfo(cool1, 'nchroms') == cool1.chroms().shape[0],\
        f'coolerInfo() failed to retrieve metadata \'nchroms\' from {fmcool1}'
    assert coolerInfo(cool2, 'nchroms') == cool2.chroms().shape[0],\
        f'coolerInfo() failed to retrieve metadata \'nchroms\' from {fmcool2}'
