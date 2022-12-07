import numpy as np
from fluctreconpy import reform_matrix, undulation_matrix

def test_undulation_matrix():
    h_matrix_ref=reform_matrix(undulation_matrix(14110, floating=0),[64, 64])
    while 1 :
        vector=(np.random.rand(64)-0.5)*4.
        scalar=np.matmul(vector, np.matmul(h_matrix_ref, vector))
        if scalar < 0 :
            print(scalar)
