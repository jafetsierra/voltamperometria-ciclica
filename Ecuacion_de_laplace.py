##---Jafet Israel Sierra---##
##---Universidad Nacional de Colombia---##
##---Resolucion numerica de la ecuacion de laplace---##
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d

P = np.zeros((100,100))

for n in range(1,40,1):
    for n in range(1,98,1):
        for m in range(1,98,1):
            if ( 35>=n>=33 and 48<=m<=50 ):
                P[n,m] = 1
            elif(68<=n<=70 and 48<=m<=50):
                P[n,m] = -1
            elif( n ==98 or m ==98 ):
                P[n,m] = 0
            else:
                P[n,m] = (P[n+1,m] + P[n-1,m] + P[n,m+1] +P[n,m-1])/4

x = np.arange(0,100,1)
y = np.arange(0,100,1)
X, Y = np.meshgrid(x,y)
z = np.array([P[x,y] for x,y in zip(np.ravel(X), np.ravel(Y))])
Z = z.reshape(X.shape)

fig = plt.figure(1, figsize=(10,10))
ax = fig.gca(projection='3d')
surf = ax.plot_surface(X,Y,Z, rstride=3,cstride=3, cmap='Blues', alpha=0.8,linewidth=0.7)
cset = ax.contourf(X,Y,Z,zdir ='x', offset=0,cmap='Oranges')
#cset = ax.contourf(X,Y,Z,zdir='y', offset=100, cmap='Blues')
cset = ax.contourf(X,Y,Z,zdir ='z', offset=-1.5, cmap='Purples')
ax.set_xlim()
ax.set_ylim()
ax.set_zlim(-1.5,1.5)
ax.set_title('Delta V en posiciones adyancentes a electrodos ')
ax.set_zlabel('Potencial (Volts)')
ax.set_xlabel('nm')
ax.set_ylabel('nm')
fig.colorbar(surf, shrink =0.5, aspect=10)
plt.show()
