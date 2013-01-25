import numpy as np
import matplotlib.pyplot as plt

font = {'size'   : 20, 'weight':'bold'}
plt.rc('xtick', labelsize=12)
plt.rc('ytick', labelsize=12)
plt.rc('font',**font)

impath = '../../resources/images/simul/'
SAVE_PLOTS = False
SHOW_PLOTS = True
n = 1e5

if n == 1e4:
    M = np.loadtxt('acccorr')
else:
    M = np.loadtxt('acccorr1e5')
M = np.transpose(M)

axx, exx = M[2], M[3]
axy, exy = M[4], M[5]
ayx, eyx = M[6], M[7]
a12, e12 = M[8], M[9]
a34, e34 = M[10], M[11]
a14, e14 = M[12], M[13]
a23, e23 = M[14], M[15]

fig1 = plt.figure()
fig1.suptitle(r'$\langle u_{x,x} u_{x,x}\rangle (= \langle u_{y,y} u_{y,y}\rangle)$')
ax1 = fig1.add_subplot(1,2,1)
ax2 = fig1.add_subplot(1,2,2)
ax1.plot(axx,'*-k',linewidth=3,markersize=10)
ax1.plot(exx,'o-r',linewidth=3)
ax2.plot(np.log10(abs(axx)),'*-k',linewidth=3,markersize=10)
ax2.plot(np.log10(abs(exx)),'o-g',linewidth=3)

fig2 = plt.figure()
fig2.suptitle(r'$\langle u_{x,y} u_{x,y}\rangle$')
ax1 = fig2.add_subplot(1,2,1)
ax2 = fig2.add_subplot(1,2,2)
ax1.plot(axy,'*-k',linewidth=3,markersize=10)
ax1.plot(exy,'o-r',linewidth=3)
ax2.plot(np.log10(abs(axy)),'*-k',linewidth=3,markersize=10)
ax2.plot(np.log10(abs(exy)),'o-g',linewidth=3)

fig3 = plt.figure()
fig3.suptitle(r'$\langle u_{y,x} u_{y,x}\rangle $')
ax1 = fig3.add_subplot(1,2,1)
ax2 = fig3.add_subplot(1,2,2)
ax1.plot(ayx,'*-k',linewidth=3,markersize=10)
ax1.plot(eyx,'o-r',linewidth=3)
ax2.plot(np.log10(abs(ayx)),'*-k',linewidth=3,markersize=10)
ax2.plot(np.log10(abs(eyx)),'o-g',linewidth=3)

fig4 = plt.figure()
fig4.suptitle(r'$\langle u_{x,x} u_{x,y}\rangle$')
ax1 = fig4.add_subplot(1,2,1)
ax2 = fig4.add_subplot(1,2,2)
ax1.plot(a12,'*-k',linewidth=3,markersize=10)
ax1.plot(e12,'o-r',linewidth=3)
ax2.plot(np.log10(abs(a12)),'*-k',linewidth=3,markersize=10)
ax2.plot(np.log10(abs(e12)),'o-g',linewidth=3)

fig5 = plt.figure()
fig5.suptitle(r'$\langle u_{y,x} u_{y,y}\rangle$')
ax1 = fig5.add_subplot(1,2,1)
ax2 = fig5.add_subplot(1,2,2)
ax1.plot(a34,'*-k',linewidth=3,markersize=10)
ax1.plot(e34,'o-r',linewidth=3)
ax2.plot(np.log10(abs(a34)),'*-k',linewidth=3,markersize=10)
ax2.plot(np.log10(abs(e34)),'o-g',linewidth=3)

fig6 = plt.figure()
fig6.suptitle(r'$\langle u_{x,x} u_{y,y}\rangle$')
ax1 = fig6.add_subplot(1,2,1)
ax2 = fig6.add_subplot(1,2,2)
ax1.plot(a14,'*-k',linewidth=3,markersize=10)
ax1.plot(e14,'o-r',linewidth=3)
ax2.plot(np.log10(abs(a14)),'*-k',linewidth=3,markersize=10)
ax2.plot(np.log10(abs(e14)),'o-g',linewidth=3)

fig7 = plt.figure()
fig7.suptitle(r'$\langle u_{x,y} u_{y,x}\rangle$')
ax1 = fig7.add_subplot(1,2,1)
ax2 = fig7.add_subplot(1,2,2)
ax1.plot(a23,'*-k',linewidth=3,markersize=10)
ax1.plot(e23,'o-r',linewidth=3)
ax2.plot(np.log10(abs(a23)),'*-k',linewidth=3,markersize=10)
ax2.plot(np.log10(abs(e23)),'o-g',linewidth=3)

if SAVE_PLOTS:
    ext = '.png'
    fig1.savefig(impath + 'Axxcorr' + str(n) + ext, format='png')
    fig2.savefig(impath + 'Axycorr' + str(n) + ext, format='png')
    fig3.savefig(impath + 'Ayxcorr' + str(n) + ext, format='png')
    fig4.savefig(impath + 'A12corr' + str(n) + ext, format='png')
    fig5.savefig(impath + 'A34corr' + str(n) + ext, format='png')
    fig6.savefig(impath + 'A14corr' + str(n) + ext, format='png')
    fig7.savefig(impath + 'A23corr' + str(n) + ext, format='png')

if SHOW_PLOTS:
    plt.show()
