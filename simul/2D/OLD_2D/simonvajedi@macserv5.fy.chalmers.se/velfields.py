import numpy as np
from utils import *
from constants import *

def turb_sfield(pos,a):
    s = 0.0
    for i in xrange( NUM_TERMS ):
        kx = i - MID
        kx *= TWO_PI
        for j in xrange( NUM_TERMS ):
            ky = j - MID
            ky *= TWO_PI
            k2 = kx*kx + ky*ky
            ang = kx*pos[0] + ky*pos[1]
            s += (a[i][j].real*np.cos( ang ) - a[i][j].imag*np.sin( ang) )* np.exp( -k2*CORRL2/4.0 )

    return SQRT_2PI * CORR_LEN * SQRT_C0 * s


def turb_vfield(pos,a):
    """
    Return the fluid velocity (tubulence vel field) for a particle 
    at position pos.
    """
    sx = 0.0
    sy = 0.0
    for i in xrange( NUM_TERMS ):
        kx = i - MID
        kx *= TWO_PI
        for j in xrange( NUM_TERMS ):
            ky = j - MID
            ky *= TWO_PI
            k2 = kx*kx + ky*ky
            ang = kx*pos[0] + ky*pos[1]
            sx += (-ky*a[i][j].real*np.sin( ang ) - ky*a[i][j].imag*np.cos( ang) )* np.exp( -k2*CORRL2/4.0 )
            sy += (-kx*a[i][j].real*np.sin( ang ) - kx*a[i][j].imag*np.cos( ang) )* np.exp( -k2*CORRL2/4.0 )

    coeff = SQRT_2PI * CORR_LEN * SQRT_C0
    return coeff * sy, -coeff*sx

def turb_acc(pos, a):
    """
    Return the matrix (vector) [u_xx, u_xy; u_yx, u_yy]
    """
    uxx, uxy, uyx, uyy = 0.0, 0.0, 0.0, 0.0
    for i in xrange( NUM_TERMS ):
        kx = i - MID
        kx *= TWO_PI
        for j in xrange( NUM_TERMS ):
            ky = j - MID
            ky *= TWO_PI
            k2 = kx*kx + ky*ky
            ang = kx*pos[0] + ky*pos[1]
            uxx += (a[i][j].real*np.cos( ang ) - a[i][j].imag*np.sin( ang) )* \
                    np.exp( -k2*CORRL2/4.0 ) * (-kx*ky)
            uxy += (a[i][j].real*np.cos( ang ) - a[i][j].imag*np.sin( ang) )* \
                    np.exp( -k2*CORRL2/4.0 ) * (-ky*ky)
            uyx += (a[i][j].real*np.cos( ang ) - a[i][j].imag*np.sin( ang) )* \
                    np.exp( -k2*CORRL2/4.0 ) * (kx*kx)
            uyy += (a[i][j].real*np.cos( ang ) - a[i][j].imag*np.sin( ang) )* \
                    np.exp( -k2*CORRL2/4.0 ) * (kx*ky)
                    
    return [uxx, uxy, uyx, uyy]

def test_turb():
    
    print STR_SEP
    print 'Test scalar corr of the turbulence field.'
    print STR_SEP

    num_avgs = 1000

    px, py = 0.0, 0.0
    
    L = 1.0
    
    vpx = np.linspace(-L/2.0,L/2.0,11)
    vpy = np.linspace(-L/2.0,L/2.0,11)

    scals = np.zeros((len(vpx),len(vpy)))
    expons = np.zeros((len(vpx),len(vpy)))

    a = [[complex(0.0,0.0) for j in xrange(NUM_TERMS)] for i in xrange(NUM_TERMS)]
    for i in xrange(len(vpx)):
        x = vpx[i]
        for j in xrange(len(vpy)):
            y = vpy[j]
            s = 0.0
            for k in xrange(num_avgs):
                compute_a(a)
                s1 = turb_sfield((px,py),a)
                s2 = turb_sfield((x,y),a)
                s += s1*s2
            s /= float(num_avgs)
            scals[i][j] = s
            expons[i][j] = np.exp( -((px-x)**2 + (py-y)**2)/2.0/CORRL2)

            print 'x =',x,'y =', y, '<phi phi> =',s, 'exp(...) =',expons[i][j]

    len2 = len(scals) * len(scals)
    np.savetxt('turbcorr',(np.reshape(scals,len2),np.reshape(scals,len2)))


