import numpy as np
import matplotlib.pyplot as plt


M = np.loadtxt('timecorr2')
M = np.transpose(M)
s = M[0]
e = M[1]
ux = M[2]
eux = M[3]
uy = M[4]
euy = M[5]

plt.figure()
plt.plot(s,'*-k',linewidth=3,markersize=10)
plt.plot(e,'^-r',linewidth=3,markersize=10)

plt.figure()
plt.plot(ux,'*-k',linewidth=3,markersize=10)
plt.plot(eux,'^-g',linewidth=3,markersize=10)

plt.figure()
plt.plot(uy,'*-k',linewidth=3,markersize=10)
plt.plot(euy,'^-b',linewidth=3,markersize=10)



plt.show()
