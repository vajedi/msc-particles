import matplotlib.pyplot as plt
import numpy as np
import random
from constants import *

def f(x):
    global TWO_PI, DSQRT_2PI, SQRT_CORRL, sigma, nTerms, mid, a
    sum = 0.0
    for i in range( nTerms ):
        k = i-mid
        k = TWO_PI * k
        sum += (a[i].real*np.cos(k*x)-a[i].imag*np.sin(k*x)) * \
                np.exp( -CORR_LENGTH2*k*k/4.0 )
    return sigma * DSQRT_2PI * SQRT_CORRL * sum
    
def set_a():
    global a, nTerms, mid
    for i in range( nTerms ):
        if i > mid:
            k = i-mid
            a[i] = a[mid-k].conjugate()
        else:
            a[i] = random.gauss(0, 1) + random.gauss(0, 1)*1j
    a[mid] = a[mid].real
    
def wrap(value, min, max):
    diff = max - min
    value = value.real
    if diff < 0:
        return value
    elif diff == 0:
        return min
    while value > max:
        value -= diff
    while value < min:
        value += diff
    return value
    
def clamp(value, min, max):
    if value < min:
        return min
    if value > max:
        return max
    return value

def compute_lya():
    global ifig, dt, nTimeSteps, t, dx0, lya_k, alpha

    nParticles = 2 #input('Nbr of particles? ')
    
    if PLOT_LYAK:
		print 'Init plots.'
		fig_parts = plt.figure( ifig )
		ifig += 1
		fig_lyak = plt.figure( ifig )
		ifig += 1

		ax_parts = fig_parts.add_subplot(1,1,1)
		ax_lyak  = fig_lyak.add_subplot(1,1,1)
		
		ax_parts.set_xlabel(r'$t$', fontsize=fontsize)
		ax_parts.set_ylabel(r'$x$', fontsize=fontsize)
		ax_parts.set_ylim([0.0,1.0])
		
		ax_lyak.set_xlabel(r'$k$', fontsize=fontsize)
		ax_lyak.set_ylabel(r'$\lambda_k$', fontsize=fontsize)
		print 'Compute particle tracks.'
    
    
    x = [[0 for i in t] for j in range(nParticles)]
    for i in range(len(x)):
        x[i][0] = initRange[0] + (initRange[1]-initRange[0])*random.random()
    
    x[1][0] = x[0][0] - dx0
        
    k = 0
    prevdx = dx0
    for i in range( len(x[0])-1 ):
        if i % NBR_OF_STEPS_PER_K == 0:
            dx = abs( x[1][i] - x[0][i] )
            if dx == 0:
                dx = 1e-10
            alpha[k] = float(dx) / prevdx
            lya_k[k] = sum( np.log(alpha[:(k+1)]) ) / (k+1)
            x[1][i] = x[0][i] - dx0
            k += 1
            prevdx = dx
        
        set_a()
        for iParticle in range(nParticles):
            ff = f( x[iParticle][i] )
            xx = x[iParticle][i]
            x[iParticle][i+1] = wrap( xx + ff, 0.0, 1.0 )

    if PLOT_LYAK and SHOW_PLOTS:
        print 'Plot the random walks.'
        for iP in range(nParticles):
            ax_parts.plot(t,x[iP],',')
        ax_lyak.plot(range(NBR_OF_K_VALUES), lya_k)
    
    nAvgs = 5 
    return sum( lya_k[-nAvgs:] ) / nAvgs

#k =     [0 for i in range(NBR_OF_K_VALUES)]
alpha = [0 for i in range(NBR_OF_K_VALUES)]
lya_k = [0 for i in range(NBR_OF_K_VALUES)]
lya   = [0 for i in range(NBR_OF_LYAK)]

ifig = 1

nTerms = 12 + 1
mid = (nTerms-1)/2
a = [0 for i in range( nTerms )]
t = [i*dt for i in range( nTimeSteps )]
var_list = np.linspace(0, 3.0*CORR_LENGTH2, NBR_OF_LYAK)#[j*dVAR for j in range(NBR_OF_LYAK)]
sigma_list = np.linspace(0, 3.0*CORR_LENGTH, NBR_OF_LYAK)
variance = VARIANCE
sigma = SIGMA

if not PLOT_LYAK:
	fig = plt.figure( ifig )
	ifig += 1
	ax = fig.add_subplot(1,1,1)

#########################################################
#
#     Compute the Lyapunov exponent.
#
#########################################################

print '\n***********************************************************\n'
print '\tCompute the Lyapunov exponent.'
print '___________________________________________________________\n\n'

for i_lya in range(NBR_OF_LYAK):
    if NBR_OF_LYAK > 1:
        #variance = var_list[i_lya]
        sigma = sigma_list[i_lya]
        print '{:.1f} %'.format( 100.0 * i_lya / NBR_OF_LYAK )
    lya[i_lya] = compute_lya()

print '\n***********************************************************\n'            

if not PLOT_LYAK:
    ax.plot(np.divide(sigma_list, CORR_LENGTH), lya)
    ax.set_ylabel(r'$\lambda$', fontsize=fontsize)
    ax.set_xlabel(r'$\sigma/\xi$', fontsize=fontsize)

if SHOW_PLOTS:
    plt.show()
