import matplotlib.pyplot as plt
import numpy as np
import random

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

def task1():
    global plt
    #########################################################
    #
    #          Task 1.
    #
    #########################################################
    print '\n***********************************************************\n'
    print '\nTask 1\n\t--Check if f(x) is Gaussian.\n'
    print '___________________________________________________________\n'


    fig = plt.figure()
    

    F = [0 for i in range(nPointsHist)]
    for i in range(len(F)):
        
        set_a()
        F[i] = f( random.random() )
    
    print 'Plot the histogram.'
    plt.hist(F, bins=60)
    print '\n***********************************************************\n'
    plt.show()

def task2():
    global plt
    #########################################################
    #
    #          Task 2.
    #
    #########################################################
    print '\n***********************************************************\n'
    print '\nTask 2\n\t--Check the correlation function of f(x).\n'
    print '___________________________________________________________\n'


    fig = plt.figure()
    

    F = [0 for i in range(nPointsHist)]
    G = [0 for i in range(nPointsHist)]
    X = [0 for i in range(nPointsHist)]
    for i in range(len(F)):
        set_a()
        x1 = random.random()
        x2 = random.random()
        f1 = f( x1 )
        f2 = f( x2 )            
        F[i] = f1 * f2
        diff = x1 - x2
        G[i] = variance*np.exp(-diff*diff/2.0/CORR_LENGTH2)
	X[i] = diff
    
    print 'Plot histograms.'
    plt.hist(F, bins=100, normed=True)
    plt.xlabel(r'$\langle f_l(x) f_{l}(x)\rangle $')
    
    fig = plt.figure()
    plt.plot(X,G, 'ok')
    plt.xlabel(r'$\sigma^2 e^{-(x-x\')^2/2\xi^2}$')
    print '\n***********************************************************\n'


def task3():
    global dt, nTimeSteps, t, ANIM_PLOTS, nUpdates
    #########################################################
    #
    #          Task 3.
    #
    #########################################################
    print '\n***********************************************************\n'
    print 'Task 3\n\t--Simulate the particles.'
    print '___________________________________________________________\n'

    nParticles = input('Nbr of particles? ')

    fig = plt.figure()

    print '\nInit starting points.'        
    
    x = [[0 for i in t] for j in range(nParticles)]
    for i in range(len(x)):
        x[i][0] = initRange[0] + (initRange[1]-initRange[0])*random.random()

    if ANIM_PLOTS:
        ax = fig.add_subplot(111)
        ax.set_xlim(0, nTimeSteps*dt)
        ax.set_ylim(0.0, 1.0)
        lines = [None for i in range(nParticles)]
        for i in range(nParticles):
            lines[i], = ax.plot(t, x[i], ',')
        
    nSkips = nTimeSteps / nUpdates
        
    print 'Compute the particle tracks.'
    for i in range( len(x[0])-1 ):
        set_a()
        for iParticle in range(nParticles):
            ff = f( x[iParticle][i] )
            xx = x[iParticle][i]
            x[iParticle][i+1] = wrap( xx + ff, 0.0, 1.0 )
        if i % nSkips == 0 and ANIM_PLOTS:
            for iParticle in range(nParticles):
                lines[iParticle].set_ydata( x[iParticle] )
            fig.canvas.draw()

    if not ANIM_PLOTS:
        print 'Plot the random walks.'
        for iP in range(nParticles):
            plt.plot(t,x[iP],',')

    print '\n***********************************************************\n'
   

def task4():
    #########################################################
    #
    #          Task 4.
    #
    #########################################################
    print '\n***********************************************************\n'
    print 'Task 4\n\t--Locate the phase transition.'
    print '___________________________________________________________\n'



#########################################################
#
#   Constants
#
#########################################################

CORR_LENGTH = 0.1
variance = 0.01
nTimeSteps = 6000
nPointsHist = 1500
initRange = [0.1, 0.9]
SHOW_PLOTS = True
ANIM_PLOTS = False
nUpdates = 60

CORR_LENGTH2 = CORR_LENGTH * CORR_LENGTH
TWO_PI = 2*np.pi
SQRT_2PI = np.sqrt( 2*np.pi )
DSQRT_2PI = np.sqrt( SQRT_2PI )
SQRT_CORRL = np.sqrt( CORR_LENGTH )
dt = 1.0
t = [i*dt for i in range( nTimeSteps )]
sigma = np.sqrt( variance )

if ANIM_PLOTS:
    plt.ion()


nTerms = 10 + 1
mid = (nTerms-1)/2
a = [0 for i in range( nTerms )]


#########################################################
#
#     Execute tasks
#
#########################################################

#task1()
#task2()
task3()


print

if SHOW_PLOTS:
    plt.show()
