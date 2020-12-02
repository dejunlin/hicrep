#! /usr/bin/env python3
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright 2020 dejun <dejun.lin@gmail.com>
# Created at 2020-12-01 19:54 on dejun@dejun-GS60-2QE
# Usage: test_sccByDiag.py
# Description: test sccByDiag
#
# Distributed under terms of the GNU General Public License v3.0.
import numpy as np
import scipy.sparse as sp
from scipy.sparse import coo_matrix
from hicrep.hicrep import (
    sccByDiag
    )

def testSccByDiag():
    m1 = coo_matrix(np.random.uniform(size=(5000, 5000)))
    m2 = coo_matrix(np.random.uniform(size=(5000, 5000)))

    scc = sccByDiag(m1, m1, m1.shape[0] - 3)

    assert np.isclose(scc, 1),\
        f"sccByDiag returns self SCC score other than one"

    scc = sccByDiag(m1, m2, m1.shape[0] - 3)
    # Test that the correlation between two random arrays is near zero
    assert scc < 1e-3 and scc > -1e-3,\
        f"sccByDiag returns SCC score significantly different than 0"\
        f"for two randomly populated matrices"
