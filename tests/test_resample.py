import numpy as np
from scipy.sparse import coo_matrix
from hicrep.utils import (
    resample
    )

def testResample():
    arr1 = coo_matrix(np.array([0, 1, 3, 1000, 3, 1, 0]))
    arr2 = coo_matrix(np.random.rand(10000) * 1000)

    #Test that resample returns the correct new number of contacts as specified by size
    assert np.sum(resample(arr1, 10000).toarray()) == 10000,\
        f"resample returns an array which doesn't sum to size (when size != 0)"
    assert np.sum(resample(arr2, 46789).toarray()) == 46789,\
        f"resample returns an array which doesn't sum to size (when size != 0)"
    assert np.sum(resample(arr1, 0).toarray()) == 0,\
        f"resample returns an array which doesn't sum to size (when size = 0)"

    #Test that the results of resample roughly match the distribution of the input
    result = resample(arr1, 100000).toarray()[0]
    #There should be no contacts in a bin which originally had no contacts
    assert result[0] == 0,\
        f"resample returns an array which doesn't sum to size (when size != 0)"
    #And there should be more contact in the bin with 1000 times more contacts
    assert result[3] > result[1],\
        f"resample returns an array which doesn't sum to size (when size = 0)"
