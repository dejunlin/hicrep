import numpy as np
from scipy.sparse import coo_matrix
from hicrep.hicrep import (
    sccOfDiag
    )

def testSccOfDiag():
    arr1 = np.random.rand(100000) * 100
    arr2 = np.random.rand(100000) * 100

    (corr, weight) = sccOfDiag(arr1, arr1)

    #Test that the correlation between a diagonal and itself is one
    assert np.isclose(corr, 1),\
        f"sccOfDiag returns correlation value other than one between a diagonal and itself"
    #Test that it returns correct weight
    assert np.isclose(weight , (100000 + 1) / 12),\
        f"sccOfDiag returns unexpected weight value"

    (corr, _) = sccOfDiag(arr1, arr2)
    #Test that the correlation between two random arrays is near zero
    assert corr < 0.25 and corr > -0.25,\
        f"sccOfDiag returns correlation significantly different than 0 for two randomly populated diagonals"
