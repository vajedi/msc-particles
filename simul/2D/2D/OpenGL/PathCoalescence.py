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

def f2D(x, y):
    global TWO_PI, variance, DSQRT_2PI, SQRT_CORRL, sigma, nTerms, mid, a
    sum = 0.0
    for i in range( nTerms ):
        kx = i-mid
        kx *= TWO_PI
        for j in range( nTerms ):
            ky = j-mid
            ky *= TWO_PI
            angle = kx*x + ky*y
            k2 = kx*kx + ky*ky
            sum += (a[i][j].real*np.cos(angle)-a[i][j].imag*np.sin(angle)) * \
                    np.exp( -CORR_LENGTH2*k2/4.0 )
    return sigma * SQRT_2PI * CORR_LENGTH * sum
    
def set_a():
    global a, nTerms, mid
    for i in range( nTerms ):
        if i > mid:
            k = i-mid
            a[i] = a[mid-k].conjugate()
        else:
            a[i] = random.gauss(0, INV_SQRT_2) + \
                   random.gauss(0, INV_SQRT_2)*1j
    a[mid] = random.gauss(0, 1)
    
def set_a2D():
    global a, nTerms, mid, INV_SQRT_2
    for i in range( nTerms ):
        for j in range( mid, nTerms ):
            re = random.gauss(0, INV_SQRT_2)
            im = random.gauss(0, INV_SQRT_2)
            a[i][j] = re + im*1j
            a[2*mid-i][2*mid-j] = re - im*1j
    a[mid][mid] = random.gauss(0, 1)
    
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
    global ifig, plt
    #########################################################
    #
    #          Task 1.
    #
    #########################################################
    print '\n***********************************************************\n'
    print '\nTask 1\n\t--Check if f(x) is Gaussian.\n'
    print '___________________________________________________________\n'


    fig = plt.figure( ifig )
    ifig += 1

    F = [0 for i in range(nPointsHist)]
    for i in range(len(F)):
        
        set_a()
        F[i] = f( random.random() )
    
    print 'Plot the histogram.'
    plt.hist(F, bins=60)
    print '\n***********************************************************\n'

def task2():
    #########################################################
    #
    #          Task 2.
    #
    #########################################################
    print '\n***********************************************************\n'
    print '\nTask 2\n\t--Check the correlation function of f(x).\n'
    print '___________________________________________________________\n'
    

    F = [0 for i in range(nPointsHist)]
    G = [0 for i in range(nPointsHist)]
    x1 = 0
    for iji in range(20):
        x2 = random.random()
        for i in range(len(F)):
            set_a()
            f1 = f( x1 )
            f2 = f( x2 )            
            F[i] = f1 * f2
        diff = x1 - x2
        G = variance*np.exp(-diff*diff/2.0/CORR_LENGTH2)

        print str(sum(F) / nPointsHist) + '     ' + str(G)

    print '\n***********************************************************\n'

def twoDim():
    global variance, a
    #########################################################
    #
    #          Check 2D.
    #
    #########################################################
    print '\n***********************************************************\n'
    print '\nCheck 2D case..\n'
    print '___________________________________________________________\n'
    
    a = [[0 for i in range( nTerms )] for j in range( nTerms )]

    F = [0 for i in range(nPointsHist)]
    G = [0 for i in range(nPointsHist)]
    
    x = 0
    n = 1000
    dx = 0.1
    s = 0.0
    while x < 1.0:
        y = -0.5
        while y < 0.5:
            for i in range(n):
                set_a2D()
                s1 = f2D(0, 0)
                s2 = f2D(x, y)
                
                s += s1*s2 / n
            expon = variance * np.exp(-(x*x+y*y)/2.0/CORR_LENGTH2)
            print '{:0}   {:1}   {:2}   {:3}'.format(x, y, s, expon)
            s = 0.0
            
            y += dx
        x += dx
    '''
    for iji in range(20):
        x1 = random.random()
        x2 = random.random()
        for i in range(len(F)):
            set_a()
            f1 = f( x1 )
            f2 = f( x2 )            
            F[i] = f1 * f2
        diff = x1 - x2
        G = variance*np.exp(-diff*diff/2.0/CORR_LENGTH2)

        print str(sum(F) / nPointsHist) + '     ' + str(G)
    '''
    print '\n***********************************************************\n'


def task3():
    global ifig, dt, nTimeSteps, t, ANIM_PLOTS, nUpdates
    #########################################################
    #
    #          Task 3.
    #
    #########################################################
    print '\n***********************************************************\n'
    print 'Task 3\n\t--Simulate the particles.'
    print '___________________________________________________________\n'

    nParticles = input('Nbr of particles? ')

    fig = plt.figure( ifig )
    ifig += 1

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
    
    
#########################################################
#
#   Constants
#
#########################################################

CORR_LENGTH = 0.1
INV_SQRT_2 = 1.0/np.sqrt(2.0)#0.7071067811865475244008443621048490392848359376884740
variance = 0.34
nTimeSteps = 6000
nPointsHist = 10000
initRange = [0.1, 0.9]
SHOW_PLOTS = True
ANIM_PLOTS = True
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



ifig = 1

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
twoDim()

print

if SHOW_PLOTS:
    plt.show()
