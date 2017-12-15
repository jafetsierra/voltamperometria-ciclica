# -*- coding: utf-8 -*-
#Simulacion de un voltamperograma ciclico
#Jafet Israel Sierra Lagos

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from matplotlib import cm
# Reaccion Ox + e- <> Red

class parameters():
    def __init__(self):
        self.v = 0.4 #Velocidad de barrido (volts/segundo)
        self.E = 0 #Potencial aplicado
        self.E0 = 0 #Potencial estandar de la Reaccion
        self.Ei = 1 #Potencial inicial
        self.Ef = -1 #Potencial final
        self.alfa = 0.5 #Parametro de simetria del proceso
        self.k0 = 0.1        #Constante estandar del proceso cinetico de tranferencia de electrones
        self.C = 1E-03       #Concentracion (Moles/metro3)
        self.F = 96485.33    #Constante de faraday (Coulombs/mol)
        self.R = 8.3144598   #Constante de los gases ideales (Joules/kelvin*mol)
        self.T = 298.15      #Temperatura (Kelvin)
        self.Dox = 1E-10     #Constante de difusion especie oxidada (metro2/segundo)
        self.Dred = 1E-10    #Constante de difusion especie reducida (metro2/segundo)
        self.deltax = 1E-06  #Intervalo de espacio (metros)
        self.A = 7.85E-07 #Area del electrodo   (metros2)
        self.kd = 0 #Constante de velocidad heterogenea del proceso directo
        self.ki = 0 #Constante de velocidad heterogenea del proceso inverso
        self.ic = 0 #Corriente catodica
        self.ia = 0 #Corriente anodica
        self.i_neta = 0
        self.n = 0 #Numero de electrones transferidos
        self.dis = 0 #Distancia maxima desde la superficie del electrodo ()
        self.n_inter = 2300 #Numero de intervalos de tiempo tomados, se toman 2300 para que los valores de DmO y DmR (adimensionales) sean menores a 0.45 y
    def parametros(self):
        self.texp = float(2*(self.Ei-self.Ef)/self.v) #Tiempo total de experimentacion ()
        print ("indique el numero de electrones transferidos")
        self.n = input()
        print ("indique la distancia maxima (micrometros)")
        self.dis = input()
        self.deltat = float(self.texp/self.n_inter) #Intervalo de tiempo  (segundos)
        self.DmO = self.Dox*self.deltat/self.deltax**2 #Si este valor adimencional es mayor a 0.45 NO CONVERGE EL METODO
        self.DmR = self.Dred*self.deltat/self.deltax**2 #Si este valor adimencional es mayor a 0.45 NO CONVERGE EL METODO
Val = parameters()
Val.parametros()

tiempo = int(Val.n_inter)
Co_matriz = np.zeros([int(Val.n_inter),int(Val.dis)]) #Matriz de puntos para especie oxidada
Cr_matriz = np.zeros([int(Val.n_inter),int(Val.dis)]) #Matriz de puntos para especie reducida
Co_i = np.zeros(int(Val.dis))
Co_f = np.zeros(int(Val.dis))
Cr_i = np.zeros(int(Val.dis))
Cr_f = np.zeros(int(Val.dis))
corriente_c = []
corriente_a = []
corriente_neta = []
potencial = []

for i in range(int(Val.dis)):
    Co_i[i] = Val.C

Co_matriz[0,:] = Co_i
Cr_matriz[0,:] = Cr_i
print (Val.n_inter)
for t in range(Val.n_inter):
    if t<Val.n_inter/2:
        Val.E = float(Val.Ei - Val.v*t*Val.deltat)
    else:
        Val.E = float(Val.Ef + Val.v*((t-tiempo/2)*Val.deltat))
    print (Val.E)
    if t>=1:
        Val.kd = float(Val.k0*np.exp((-Val.alfa*Val.F*Val.n*(Val.E-Val.E0))/(Val.R*Val.T))) #Calculo de la constante cinetica directa
        Val.ki = float(Val.k0*np.exp(((1-Val.alfa)*Val.F*Val.n*(Val.E-Val.E0))/(Val.R*Val.T))) #Calculo de la constante cinetica inversa
        for x in range(Val.dis-1):
            if x==0:
                Co_f[x] = (Co_i[x] + (-Val.kd*Co_i[x] + Val.ki*Cr_i[x])*Val.deltat) #Calculo de [] especie oxidada por el proceso cinetico de transferencia de electrones
                Cr_f[x] = (Cr_i[x] + (Val.kd*Co_i[x] - Val.ki*Cr_i[x])*Val.deltat) #Calculo de [] especie reducida por el proceso cinetico de transferencia de electrones
                if Co_f[x] > Val.C:
                    Co_f[x] = Val.C
                if Co_f[x] <0:
                    Co_f[x] =0
                if Cr_f[x] <0:
                    Cr_f[x] = 0
            else:
                Co_f[x] = (Co_i[x] + Val.DmO*(Co_i[x+1]-2*Co_i[x]+Co_i[x-1])) #discretizacion segunda ley de fick
                Cr_f[x] = (Cr_i[x] + Val.DmR*(Cr_i[x+1]-2*Cr_i[x]+Cr_i[x-1]))  #discretizacion segunda ley de fick
        Co_f[Val.dis-1] = Val.C #Condicion de frontera a distancia "infinita"
        Cr_f[Val.dis-1] = 0 #Condicion de frontera a distancia "infinita"
        Co_matriz[t,:] = Co_f
        Cr_matriz[t,:] = Cr_f
        #Calculo corriente generada por el transporte de masa (difusion)
        Val.ic = (Val.n *Val.F*Val.A*Val.Dox*((Co_f[2]-Co_f[0])*1E09/Val.deltax)) #Nanoampere
        Val.ia = (-Val.n *Val.F*Val.A*Val.Dred*((Cr_f[2]-Cr_f[0])*1E09/Val.deltax)) #Nanoampere
        potencial.append(float(Val.E))
        corriente_c.append(float(Val.ic))
        corriente_a.append(float(Val.ia))
        Co_i = Co_f
        Cr_i = Cr_f
for i in range(len(corriente_a)):
    corriente_neta.append(corriente_a[i]+corriente_c[i])
print (len(corriente_neta))
####################################Grafica para la expecie oxidada
fig1 = plt.figure(1)
ax1 = fig1.gca(projection='3d')
x = np.arange(1,Val.n_inter,1)
y = np.arange(1,Val.dis-1,1)
X,Y = np.meshgrid(x,y)
z1 = np.array([Co_matriz[x,y] for x,y in zip(np.ravel(X), np.ravel(Y))])
Z1 = z1.reshape(X.shape)
c1 = int(Val.n_inter/100)
c2 = int(Val.dis/10)

surf = ax1.plot_surface(X,Y,Z1, rstride=c1, cstride=c2, alpha='0.8',cmap=cm.coolwarm,linewidth=0.9)
cset = ax1.contour(X,Y,Z1,zdir='x',offset=0,cmap='Reds')
cset = ax1.contour(X,Y,Z1,zdir='y',offset=Val.dis,cmap='Reds')
cset = ax1.contour(X,Y,Z1,zdir='z',offset=0,cmap='Purples')
ax1.set_title('[] especie oxidada como funcion de la posicion y el tiempo ')
ax1.set_xlabel('unidad de Tiempo (cada unidad es de 4.35 ms)')
ax1.set_ylabel('Distancia (um)')
ax1.set_zlabel('Concentracion (mol/m^3)')
ax1.set_xlim(0,Val.n_inter)
ax1.set_ylim(0,Val.dis)
ax1.set_zlim(0,Val.C)
fig1.colorbar(surf, shrink =0.5, aspect=10)
###################################Grafica para la expecie reducida
fig2 = plt.figure(2)
ax2 = fig2.gca(projection='3d')
z2 = np.array([Cr_matriz[x,y] for x,y in zip(np.ravel(X), np.ravel(Y))])
Z2 = z2.reshape(X.shape)

surf = ax2.plot_surface(X,Y,Z2, rstride=c1, cstride=c2, alpha='0.8',cmap=cm.coolwarm,linewidth=0.9)
cset = ax2.contour(X,Y,Z2,zdir='x',offset=0,cmap='Reds')
cset = ax2.contour(X,Y,Z2,zdir='y',offset=Val.dis,cmap='Reds')
cset = ax2.contour(X,Y,Z2,zdir='z',offset=0,cmap='Purples')
ax2.set_title('[] especie reducida como funcion de la posicion y el tiempo ')
ax2.set_xlabel('unidad de Tiempo (cada unidad es de 4.35 ms)')
ax2.set_ylabel('Distancia (um)')
ax2.set_zlabel('Concentracion (mol/m^3)')
ax2.set_xlim(0,Val.n_inter)
ax2.set_ylim(0,Val.dis)
ax2.set_zlim(0,Val.C)
fig2.colorbar(surf, shrink =0.5, aspect=10)
###################################Grafica de corriente vs potencial aplicado
fig3 = plt.figure(3)
ax3 = fig3.gca()
ax3.plot(potencial,corriente_neta)
ax3.set_xlabel('Potencial (V)')
ax3.set_ylabel('Corriente (nA)')
ax3.set_title('Voltamperograma Ciclico')
ax3.grid(True)
plt.show()
###################################Paso opcional creacion de archivos de texto con los datos
import pandas as pd
libro = pd.DataFrame(data=corriente_a, columns=['corriente'])
libro.to_csv('datos_corriente.text', index=False, header=False, sep=',')

libro = pd.DataFrame(data=potencial, columns=['corriente'])
libro.to_csv('datos_potencial.text', index=False, header=False, sep=',')
