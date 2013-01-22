import numpy as np

CORR_LENGTH = 0.3
VARIANCE = 1e-6#0.01*CORR_LENGTH**2#0.00003
SIGMA = np.sqrt( VARIANCE )

VERBOSE = False

################
#  Lyapunov
################
dx0 = 0.02
NBR_OF_K_VALUES = 200
NBR_OF_STEPS_PER_K = 200
VAR_MAX = 8*CORR_LENGTH*CORR_LENGTH
NBR_OF_LYAK = 20
PLOT_LYAK = False
if (NBR_OF_LYAK == 1):
	PLOT_LYAK = True
dVAR = VARIANCE / 10.0
#################

nTimeSteps = NBR_OF_STEPS_PER_K * NBR_OF_K_VALUES

TWO_PI = 2*np.pi
CORR_LENGTH2 = CORR_LENGTH * CORR_LENGTH
SQRT_2PI = np.sqrt( 2*np.pi )
DSQRT_2PI = np.sqrt( SQRT_2PI )
SQRT_CORRL = np.sqrt( CORR_LENGTH )
dt = 1.0

initRange = [0.4, 0.6]

SHOW_PLOTS = True
fontsize = 20
