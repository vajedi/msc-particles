import numpy as np
import random
from constants import *

def compute_a(a):
    for i in xrange( NUM_TERMS ):
        for j in xrange( MID,NUM_TERMS ):
            #a[i][j] = random.gauss(0.0,INV_SQRT_2) + 1j*random.gauss(0.0,INV_SQRT_2)
            a[i][j] = a[i][j] * (1 - TIME_RATE) + \
                      complex(random.gauss(0.0,std_dev), 
                              random.gauss(0.0,std_dev))
            a[2*MID-i][2*MID-j] = complex(a[i][j].real,-a[i][j].imag)
    a[MID][MID]= complex(a[MID][MID] * (1 - TIME_RATE) + 
                  random.gauss(0.0, 1.0), 0)           
    #random.gauss(0.0, TIME_DEV_SQRT2),0)

    #return a
