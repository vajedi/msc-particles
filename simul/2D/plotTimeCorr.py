import numpy as np
import matplotlib.pyplot as plt


M = np.loadtxt('timecorr')
M = np.transpose(M)
tr = M[0]
ti = M[1]
e = M[2]

plt.figure()
plt.plot(tr,'*-k',linewidth=3,markersize=10)
plt.plot(ti,'^-g',linewidth=3,markersize=10)
plt.plot(e,'o-r',linewidth=3)


plt.show()
