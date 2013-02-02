import numpy as np
import matplotlib.pyplot as plt
import random

def randgauss():
    S = 0
    #print U, V
    while S == 0 or S >= 1:
        U = 2*np.random.random()-1
        V = 2*np.random.random()-1
        S = U*U + V*V
    return U*np.sqrt(-2.0*np.log(S)/S)

def randg():
    return np.sqrt(-2.0*np.log(random.random()))*np.cos(2*np.pi*random.random())

nPts = 10000
nTerms = 121

vec = [randgauss() for i in range(nPts)]

print sum([sum([random.gauss(0,1) for i in range(10)])/nPts for j in range(nPts)])

'''
fig = plt.figure(1)

plt.hist(vec, bins=100, normed=True)

fig = plt.figure(2)

plt.hist([random.gauss(0,1) for i in range(nPts)], bins=100, normed=True)

plt.show()
'''
