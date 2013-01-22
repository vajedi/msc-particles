import numpy as np
import random

INV_SQRT_2 = 1.0 / np.sqrt(2)


def compute_a(a, num_terms, mid):
    """
    Compute the rank-2 tensor a. a is a 
    num_terms x num_terms matrix.


    a         -- the destination matrix
    num_terms -- the number of terms to take into account 
                  when computing a.
    """

    mid = num_terms / 2

    for i in xrange(num_terms):
        for j in xrange(mid, num_terms):
            a[i][j] = random.gaussian(0,1) + 1j*
                      random.gaussian(0,1)

            a[2*mid-i][2*mid-j] = a[i][j].conjugate

    a[mid][mid] = random.gaussian(0,1)
           

def turb_field(x, a, num_terms, mid, sigma):
    """
    The scalar field f that simulate turbulence 
    (using matrix a) for a particle at position x. 


    """


