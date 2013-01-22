import numpy as np
import matplotlib.pyplot as plt


M = np.loadtxt('scalarcorr')

M = np.transpose(M)

pp = M[2]
ee = M[3]


plt.plot(np.log(np.abs(pp)),'*-k',linewidth=3,markersize=10)
plt.plot(np.log(np.abs(ee)),'o-r',linewidth=3)

plt.show()
