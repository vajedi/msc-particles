import numpy as np

NUM_TERMS = 20+1
MID = (NUM_TERMS - 1) / 2

CORR_LEN = 0.1
CORRL2 = CORR_LEN * CORR_LEN

dt = 0.001
CORR_TIME = 1.0
dt = CORR_TIME 
TIME_RATE = dt / CORR_TIME
TIME_DEV = np.sqrt( dt / CORR_TIME )
TIME_DEV_SQRT2 = np.sqrt( 2 * dt / CORR_TIME )

INV_SQRT_2 = 1.0 / np.sqrt( 2 )

std_dev = INV_SQRT_2 * TIME_DEV_SQRT2
std_dev = INV_SQRT_2

STD_DEV = 1.0

C0 = 1
SQRT_C0 = np.sqrt( C0 )

TWO_PI = 2*np.pi
SQRT_2PI = np.sqrt( TWO_PI )


STR_SEP = '\n********************************************************\n'
