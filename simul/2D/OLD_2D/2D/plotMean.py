
from numpy import pi, sin, cos, mgrid
import numpy as np

def transpose(mat):
    newmat = np.zeros((len(mat[0]), len(mat)))
    for i in range(len(newmat)):
        for j in range(len(newmat[0])):
            newmat[i][j] = mat[j][i]
    return newmat

sc = np.loadtxt( 'scalarcorr' )
x = [i[0] for i in sc]
y = [i[1] for i in sc]
scf = [i[2] for i in sc]
e = [i[3] for i in sc]

lst = [-0.5 + i*0.1 for i in range(10)]
X, Y = np.meshgrid(lst, lst)

S = np.zeros((len(lst),len(lst)))
E = np.zeros((len(lst),len(lst)))
for i in range(len(lst)):
    for j in range(len(lst)):
        S[j][i] = scf[len(lst)*i + j]
        E[i][j] = e[len(lst)*i + j]
#print x, y, scf, e

#X, Y = np.meshgrid(x, y)

X = transpose(X)
Y = transpose(Y)
S = transpose(S)

from mayavi import mlab
#s = mlab.mesh(x, y, scf)
mlab.figure(bgcolor=(1,1,1))
s = mlab.surf(X, Y, np.subtract(S,E), warp_scale="auto")
mlab.axes(nb_labels=5, line_width=2., xlabel="x", ylabel="y", zlabel="z", x_axis_visibility=True)
#t = mlab.surf(X, Y, E)
mlab.show()
