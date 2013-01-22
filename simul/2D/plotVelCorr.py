import numpy as np
import matplotlib.pyplot as plt


M = np.loadtxt('velcorr')

M = np.transpose(M)

ux = M[2]
ex = M[4]

uy = M[3]
ey = M[5]

plt.figure()
plt.plot(np.log(np.abs(ux)),'*-k',linewidth=3,markersize=10)
plt.plot(np.log(np.abs(ex)),'o-r',linewidth=3)

plt.figure()
plt.plot(np.log(np.abs(uy)),'*-k',linewidth=3,markersize=10)
plt.plot(np.log(np.abs(ey)),'o-g',linewidth=3)


plt.show()
