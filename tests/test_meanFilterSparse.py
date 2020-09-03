from scipy.sparse import coo_matrix
import numpy as np
from hicrep.utils import (
    meanFilterSparse
    )

def testMeanFilterSparse():
    #[[ 1  2  3  4]
    #[ 5  6  7  8]
    #[ 9 10 11 12]
    #[13 14 15 16]]
    mat1 = np.arange(1, 17).reshape(4,4)

    #Test that smaller values surrounded by larger ones are increased by smoothing
    assert meanFilterSparse(coo_matrix(mat1), 1).toarray()[0][0] > mat1[0][0],\
        f"Smoothing by meanFilterSparse doesn't increase relatively small values"

    #Test that larger values surrounded by smaller ones are decreased by smoothing
    assert meanFilterSparse(coo_matrix(mat1), 1).toarray()[3][3] < mat1[3][3],\
        f"Smoothing by meanFilterSparse doesn't decrease relatively large values"

    #Test that values surrounded by equally smaller and larger values are unchanged
    assert meanFilterSparse(coo_matrix(mat1), 1).toarray()[1][1] == mat1[1][1],\
        f"Smoothing by meanFilterSparsechanges values that shouldn't be changed"

    mat2 = np.matrix([[ 0,  0,  0,  0,  0],
                      [ 0,  0,  0,  0,  0],
                      [ 0,  0,  9,  0,  0],
                      [ 0,  0,  0,  0,  0],
                      [ 0,  0,  0,  0,  0]])

    mat2_smoothed = np.matrix([[ 0,  0,  0,  0,  0],
                      [ 0,  1,  1,  1,  0],
                      [ 0,  1,  1,  1,  0],
                      [ 0,  1,  1,  1,  0],
                      [ 0,  0,  0,  0,  0]])

    #Test that smoothing a matrix with a single 9 in its center using a window size of 3
    #gives a matrix with nine ones in the bins around the center
    assert np.array_equal(meanFilterSparse(coo_matrix(mat2), 1).toarray(), mat2_smoothed),\
        f"Smoothing on matrix with single non-zero entry gives unexpected result"



testMeanFilterSparse()
