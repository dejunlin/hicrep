#! /usr/bin/env python3
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright 2020 dejun <dejun.lin@gmail.com>
# Created at 2020-12-01 22:25 on dejun@dejun-GS60-2QE
# Usage: test_upperDiagCsr.py
# Description: Test upperDiagCsr
#
# Distributed under terms of the GNU General Public License v3.0.
import numpy as np
import scipy.sparse as sp
from scipy.sparse import coo_matrix, dia_matrix
from hicrep.hicrep import (
    upperDiagCsr
    )

def testUpperDiagCsr():
    m = sp.random(100, 100)
    result = upperDiagCsr(m, m.shape[0] - 2).toarray()
    for i in range(result.shape[0]):
        iResult = result[i]
        iDiag = m.diagonal(i + 1)
        nPads = iResult.size - iDiag.size
        iExpected = np.concatenate((iDiag, np.zeros(nPads)))
        assert (iResult == iExpected).all(),\
            f"upperDiagCsr produces wrong result for the {i+1} diagonal"
