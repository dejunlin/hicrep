#! /usr/bin/env python3
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright 2020 dejun <dejun.lin@gmail.com>
# Created at 2020-12-01 19:02 on dejun@dejun-GS60-2QE
# Usage: test_varVstran.py
# Description: Test varVstran
#
# Distributed under terms of the GNU General Public License v3.0.
import warnings
import numpy as np
from hicrep.utils import (
    varVstran
    )

def testVarVsTran():
    ns = np.arange(1000)
    results = varVstran(ns)
    for n in ns:
        result = results[n]
        with warnings.catch_warnings(), np.errstate(divide='ignore', invalid='ignore'):
            warnings.simplefilter("ignore", category=RuntimeWarning)
            expected = np.var(np.arange(1, n + 1) / n, ddof=1)
            assert (np.isnan(result) and np.isnan(expected)) or\
                np.allclose(result, expected),\
                f"varVstran({n}) gives wrong result: {result} while "\
                f"{expected} is expected"
