# -*- coding: utf-8 -*-
"""
Created on Tue Mar 15 17:17:49 2016

@author: Admin
"""

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
import numpy as np

numero_de_tiempos=50
numero_de_distancias=24
Do_deltat_deltax=0.22
Eo=0
Delta_de_potencial=0.6
Ef=(Eo-Delta_de_potencial)
Ei=(Eo+Delta_de_potencial)
Vd=(Ei-Ef)/(numero_de_tiempos/2)
Vi=(Ef-Ei)/(numero_de_tiempos/2)
n=1
F=96500
R=8.314
T=20+273.15
nFRT=(n*F)/(R*T)
alfa=0.5
kf=10000000
mitad_de_tiempo=int((numero_de_tiempos/2))
delta_temp=0.001
# print(Vd)
# print(Vi)
# print(Ef)
# print(Ei)
# print(nFRT)
matriz_2=np.zeros((numero_de_tiempos,3))

for u in range (0,mitad_de_tiempo):
    matriz_2[u][0]=Ei-(Vd*u)
for u2 in range (mitad_de_tiempo,numero_de_tiempos):
    matriz_2[u2][0]=Ef-(Vi*(u2-(mitad_de_tiempo-1)))
for v in range (0,numero_de_tiempos):
    matriz_2[v][1]=matriz_2[v][0]-Eo
# print(matriz_2)

matriz_1=np.zeros((numero_de_tiempos,numero_de_distancias))

for x in range (0,numero_de_distancias):
    matriz_1[0][x]=1
for c in range (0,numero_de_tiempos):
        matriz_1[c][numero_de_distancias-1]=1
# print(matriz_1)


for a in range (1,numero_de_tiempos):
    for z in range (1,numero_de_distancias-1):

        matriz_1[a][0]=(matriz_1[a][1]+(kf*np.exp((1-alfa)*nFRT*matriz_2[a][1])))/((kf*((np.exp(-alfa*nFRT*matriz_2[a][1]))+(np.exp((1-alfa)*nFRT*matriz_2[a][1]))))+1)
        matriz_1[a][z]=(matriz_1[a-1][z])+Do_deltat_deltax*(matriz_1[a-1][z-1]-2*(matriz_1[a-1][z])+matriz_1[a-1][z+1])
           
for w in range (0,numero_de_tiempos):
    matriz_2[w][2]=matriz_1[w][1]-matriz_1[w][0]

# print(matriz_2)
# print(matriz_1)
intensidad=np.zeros((numero_de_tiempos,1))
potencial=np.zeros((numero_de_tiempos))
for g in range(0,numero_de_tiempos):
    intensidad[g]=matriz_2[g][2]
    potencial[g]=matriz_2[g][0]

# print(intensidad)
# print(potencial)   

grafica = plt.figure(1)
superficie = grafica.add_subplot(121, projection='3d') #INSTRUCCIONES PARA HACER UNA GRAFICA.
y = range(numero_de_tiempos)
x = range(numero_de_distancias)
X, Y = np.meshgrid(x, y) #CUADRICULA PARA GRAFICAR.
o= 0;
for i in (X):
	print (o,  i)
	o +=1


plt.ylabel('Tiempo')
plt.xlabel('Distancia')
superficie.plot_surface(X, Y, matriz_1) #LOS VALORES X, Y Y Z A GRAFICAR.


plt.figure(2)
plt.subplot(111)
plt.plot(potencial,intensidad)
plt.xlabel('intensidad')
plt.ylabel('potencial')
plt.show()