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
from scipy.sparse import coo_matrix
import numpy as np
from hicrep.utils import (
    trimDiags
    )

def testTrimDiags():
    #[[ 1  2  3  4]
    #[ 5  6  7  8]
    #[ 9 10 11 12]
    #[13 14 15 16]]
    mat1 = coo_matrix(np.arange(1, 17).reshape(4,4))

    # Testing that trimDiags does nothing when iDiagMax is greater than n and
    # keepMain is true
    assert np.array_equal(trimDiags(mat1, 10, True).toarray() , mat1.toarray()),\
        f"trimDiags failed to modify the test matrix as expected whem no"\
        f"change needed to be made"

    # Testing that trimDiags removes the center diagonal and nothing else when
    # iDiagMax is greater than n and keepMain is false
    mat2 = np.array([[ 0,  2,  3,  4],
                     [ 5,  0,  7,  8],
                     [ 9, 10,  0, 12],
                     [13, 14, 15,  0]])
    assert np.array_equal(trimDiags(mat1, 4, False).toarray() , mat2),\
        f"trimDiags failed to modify the test matrix as expected when "\
        f"removing main diagonal"

    # Testing that diagonals greater than iDiagMax are removed
    mat3 = np.array([[ 1,  2,  0,  0],
                     [ 5,  6,  7,  0],
                     [ 0, 10, 11, 12],
                     [ 0,  0, 15, 16]])
    assert np.array_equal(trimDiags(mat1, 2, True).toarray() , mat3),\
        f"trimDiags failed to modify the test matrix as expected trimming "\
        f"outer 2 diagonals"

    # Testing the same as above except with keepMain set to false
    mat4 = np.array([[ 0,  2,  0,  0],
                     [ 5,  0,  7,  0],
                     [ 0, 10,  0, 12],
                     [ 0,  0, 15,  0]])
    assert np.array_equal(trimDiags(mat1, 2, False).toarray() , mat4),\
        f"trimDiags failed to modify the test matrix as expected when "\
        f"removing main diagonal and outer two diagonals"

    # Testing that if 0 diagonals are kept, you get a matrix of zeros with or
    # without keepMain
    mat5 = np.zeros((4,4))
    assert np.array_equal(trimDiags(mat1, 0, True).toarray() , mat5),\
        f"trimDiags failed to modify the test matrix as expected when "\
        f"keeping no diagonals (with keepMain)"
    assert np.array_equal(trimDiags(mat1, 0, False).toarray() , mat5),\
        f"trimDiags failed to modify the test matrix as expected when "\
        f"keeping no diagonals (without keepMain)"
