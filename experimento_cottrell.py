#Simulacion del experimento de Cottrell
#Jafet Israel Sierra Lagos

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from matplotlib import cm
# Reaccion unicamente de reduccion, debido a la
# aplicacion de un potencial lo suficientemente grande
# (velocidad de transferencia de electrones infinita)

class constantes():
    def __init__(self):
        self.C = 1E-03       #Concentracion (Moles/metro3)
        self.F = 96485.33    #Constante de faraday (Coulombs/mol)
        self.R = 8.3144598   #Constante d elos gases ideales (Joules/kelvin*mol)
        self.T = 298.15      #Temperatura (Kelvin)
        self.Dox = 1E-10     #Constante de difusion (metro2/segundo)
        self.deltax = 1E-06  #Intervalo de espacio (metros)
        self.deltat = 1E-03  #Intervalo de tiempo  (segundos)
        self.A = 7.85E-07 #Area del electrodo   (metros2)
        self.Dm = self.Dox*self.deltat/self.deltax**2 #Si este valor adimencional es mayor a 0.45 NO CONVERGE EL METODO
        self.n = 0
        self.texp = 0
        self.dis = 0
    def parametros(self):
        print "indique el numero de electrones transferidos"
        self.n = input()
        print "indique el tiempo del experimento (milisegundos)"
        self.texp = input()
        print "indique la distancia maxima (micrometros)"
        self.dis = input()

Cons = constantes()
Cons.parametros()
C_matriz = np.zeros((Cons.texp,Cons.dis))
C_i = np.zeros(Cons.dis)
C_f = np.zeros(Cons.dis)

for i in range(Cons.dis):
    C_i[i] = Cons.C
C_matriz[0,:] = C_i
for t in range(Cons.texp):
    if t>=1:
        for x in range(Cons.dis-1):
            if x>=1 and x<=(Cons.dis):
                C_f[x] = float(C_i[x] + Cons.Dm*(C_i[x+1]-2*C_i[x]+C_i[x-1]))
            else:
                C_f[x] = 0
            C_f[Cons.dis-1] = Cons.C
        C_matriz[t,:] = C_f
        C_i = C_f
fig = plt.figure(1)
ax = fig.gca(projection='3d')
x = np.arange(1,Cons.texp,1)
y = np.arange(1,Cons.dis-1,1)
X,Y = np.meshgrid(x,y)
z = np.array([C_matriz[x,y] for x,y in zip(np.ravel(X), np.ravel(Y))])
Z = z.reshape(X.shape)
c1 = int(Cons.texp/100)
c2 = int(Cons.dis/10)
surf = ax.plot_surface(X,Y,Z, rstride=c1, cstride=c2, alpha='0.8',cmap=cm.coolwarm,linewidth=0.9)
cset = ax.contour(X,Y,Z,zdir='x',offset=0,cmap='Reds')
cset = ax.contour(X,Y,Z,zdir='y',offset=Cons.dis,cmap='Reds')
cset = ax.contour(X,Y,Z,zdir='z',offset=0,cmap='Purples')
ax.set_title('Concentracion de la especie oxidada como funcion de la posicion y el tiempo ')
ax.set_xlabel('Tiempo (ms)')
ax.set_ylabel('Distancia (um)')
ax.set_zlabel('Concentracion (mol/m^3)')
ax.set_xlim(0,Cons.texp)
ax.set_ylim(0,Cons.dis)
ax.set_zlim(0,Cons.C)
fig.colorbar(surf, shrink =0.5, aspect=10)
plt.show()
