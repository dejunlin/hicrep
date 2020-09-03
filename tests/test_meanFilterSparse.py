import scipy.sparse as sp
import scipy.ndimage as spi
from scipy.signal import convolve2d
import numpy as np
from hicrep.utils import (
    meanFilterSparse
    )

def testMeanFilterSparse():
    # test for a range of matrix sizes and filter sizes. Can't afford testing
    # size > 200 due to exponential scaling of dense matrix manipulation
    sizes = [10, 20, 50, 100, 200]
    for size in sizes:
        a =  sp.random(size, size, density=0.5)
        for h in range(1, int((size - 1) / 2)+1, int(size / 10)):
            aResults = meanFilterSparse(a, h)
            fSize = 2 * h + 1
            kSumNeighbors = np.ones((fSize, fSize), dtype=float)
            sumNeigbhors = spi.convolve(a.todense(), kSumNeighbors, mode='constant', cval=0.0)
            nNeighbors = spi.convolve(np.ones_like(sumNeigbhors, dtype=float),
                                      kSumNeighbors, mode='constant', cval=0.0)
            aExpected = sumNeigbhors / nNeighbors
            diff = aResults - sp.coo_matrix(aExpected)
            assert np.isclose(diff.data, np.zeros(diff.nnz, dtype=float),
                              rtol=1e-12, atol=1e-12).all(),\
                f"Wrong mean filter results for sparse {size}x{size} matrix:\n"\
                f"{a.toarray()}\n"\
                f"with filter size {h}. Expected answer is:\n"\
                f"{aExpected}\nBut got answer:\n"\
                f"{aResults.toarray()}\n"\
                f"Difference is:\n"\
                f"{diff.toarray()}"\
                f"Max difference is:\n"\
                f"{diff.max()}"
