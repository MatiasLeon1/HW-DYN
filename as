# -*- coding: utf-8 -*-
"""
Created on Sat May 29 10:31:17 2021

@author: Matias
"""

import matplotlib.pyplot as plt

from matplotlib.pylab import *
from matplotlib import pyplot
import numpy as np
import math

#Interpretacion de datos ------------------------------------------------------
files = ["Los Gatos.txt","Tarzana.txt"]

G = []
T = []

G2 = []
T2 = []

L = [G,T]
L2 = [G2,T2]

con=0
for i in files:
    k=L[con]
    j=L2[con]
    x=open(i)
    for i in x:
        line = i
        ce = line.split(" ")
        for b in ce:
            if b != "":
                k.append(float(b)/980)
                j.append(float(b)/100)
                
    con+=1           
     
dt = (1/50)

tg=[];tt=[]
lts = [tg,tt]

con=0
for i in lts:
    ti=0
    while len(i) != len(L[con]):
        i.append(ti)
        ti+=dt
        
    con+=1

#/// [1] //////////////////////////////////////////////////////////////////////
#PREGUNTA 1.1
#Graficar Registros de Tarzana y Los Gatos    

plt.figure()   
plt.title("Earthquake Motion Record")
plt.plot(tt, T,linewidth = 1.5,label="Tarzana 1994",alpha=0.7)
plt.plot(tg, G,linewidth = 1.5,label="Los Gatos 1989",alpha=0.7)
plt.legend(loc="upper right",fontsize="larger")
plt.grid(b=True, alpha=0.7,linestyle="dashed", color="dimgray",linewidth=0.5)
plt.xlabel("Time [s]",size = 12)
plt.ylabel("Acceleration [g]",size = 12)
plt.xlim(0,30)
plt.axhline(0, linewidth=1, color ="k")
savefig("Earthquake_Record.png")
show()

#PREGUNTA 1.2
#PREGUNTA 1.2.1
#--------------------- CHI 3% ---------------------

T12=[]
t12=0.04

while t12 <=3:
    T12.append(t12)
    t12+=0.003
    
#Newmark’s method: linear acceleration.

beta=1/6
gamma=1/2

"Condiciones Iniciales"

"2.0""Calculos para pasos de tiempo i=0,1,2,....n"

lista_u_max=[]
lista_u1p_max=[]
lista_u2p_max=[]

lista_u1p_max_pseudo=[]
lista_u2p_max_pseudo=[]

# LOS GATOS
for j in T12: 
    chi=0.03
    masa = 1.0 #kg
    w=(2*np.pi)/j
    c= 2*chi*w #kg/s 
    k=w*w*masa #kg/s2
    dt1=1/50
    
    # Condiciones iniciales
    
    u0=0
    u1p0=0
    u2p0= G2[0]/-masa

    lista_u=[u0]
    lista_u1p=[u1p0]
    lista_u2p=[u2p0]
    lista_p_g=[]
    
    a1=(1/(beta*(dt1*dt1))*masa)+(gamma/(beta*dt1))*c
    a2=(1/(beta*dt1))*masa + (gamma/beta-1)*c
    a3=((1/(2*beta))-1)*masa + dt1*(gamma/(2*beta)-1)*c

    k_g= k + a1

    for i in range(len(G2)-1):
     
        p_g=G2[i+1]+a1*lista_u[i]+a2*lista_u1p[i]+a3*lista_u2p[i]
        lista_p_g.append(p_g)
        
        u_i=lista_p_g[i]/k_g
        
        u1p_i=(gamma/(beta*dt1))*(u_i-lista_u[i])+(1-(gamma/beta))*lista_u1p[i]+dt1*(1-(gamma/(2*beta)))*lista_u2p[i]

        u2p_i=(1/(beta*(dt1*dt1)))*(u_i-lista_u[i])-(1/(beta*dt1))*lista_u1p[i]-(1/(2*beta)-1)*lista_u2p[i]

        lista_u.append(u_i)
        lista_u1p.append(u1p_i)
        lista_u2p.append(u2p_i)
        
    lista_u_max.append(abs(max(lista_u,key=abs)))
    lista_u1p_max.append(abs(max(lista_u1p,key=abs)))   
    lista_u2p_max.append(abs(max(lista_u2p,key=abs)))
    
    lista_u1p_max_pseudo.append(w*abs(max(lista_u,key=abs)))
    lista_u2p_max_pseudo.append(w*w*abs(max(lista_u,key=abs))) 

# TARZANA

listaT_u_max=[]
listaT_u1p_max=[]
listaT_u2p_max=[]

listaT_u1p_max_pseudo=[]
listaT_u2p_max_pseudo=[]

for j in T12: 
    chi=0.03
    masa = 1.0 #kg
    w=(2*np.pi)/j
    c= 2*chi*w #kg/s 
    k=w*w*masa #kg/s2
    dt1=1/50
    
    # Condiciones iniciales
    
    u0=0
    u1p0=0
    u2p0= T2[0]/-masa

    lista_u=[u0]
    lista_u1p=[u1p0]
    lista_u2p=[u2p0]
    lista_p_g=[]
    
    a1=(1/(beta*(dt1*dt1))*masa)+(gamma/(beta*dt1))*c
    a2=(1/(beta*dt1))*masa + (gamma/beta-1)*c
    a3=((1/(2*beta))-1)*masa + dt1*(gamma/(2*beta)-1)*c

    k_g= k + a1

    for i in range(len(T2)-1):
     
        p_g=T2[i+1]+a1*lista_u[i]+a2*lista_u1p[i]+a3*lista_u2p[i]
        lista_p_g.append(p_g)
        
        u_i=lista_p_g[i]/k_g
        
        u1p_i=(gamma/(beta*dt1))*(u_i-lista_u[i])+(1-(gamma/beta))*lista_u1p[i]+dt1*(1-(gamma/(2*beta)))*lista_u2p[i]

        u2p_i=(1/(beta*(dt1*dt1)))*(u_i-lista_u[i])-(1/(beta*dt1))*lista_u1p[i]-(1/(2*beta)-1)*lista_u2p[i]

        lista_u.append(u_i)
        lista_u1p.append(u1p_i)
        lista_u2p.append(u2p_i)
        
    listaT_u_max.append(abs(max(lista_u,key=abs)))
    listaT_u1p_max.append(abs(max(lista_u1p,key=abs)))
    listaT_u2p_max.append(abs(max(lista_u2p,key=abs)))    
    
    listaT_u1p_max_pseudo.append(w*abs(max(lista_u,key=abs)))
    listaT_u2p_max_pseudo.append(w*w*abs(max(lista_u,key=abs))) 
    
    
plt.figure()
plt.title("Relative Displacement Spectra")
plt.plot(T12,lista_u_max ,linewidth = 1.5,label="Los Gatos 3%",alpha=0.7,color="b")
plt.plot(T12,listaT_u_max ,linewidth = 1.5,label="Tarzana 3%",alpha=0.7,color="r")
plt.legend(loc="upper right",fontsize="medium")
plt.grid(b=True, alpha=0.7,linestyle="dashed", color="dimgray",linewidth=0.5)
plt.xlabel("Period [s]",size = 12)
plt.ylabel("Sd [m]",size = 12)
plt.xlim(0,3)
plt.ylim(0,0.9)
# plt.axhline(0, linewidth=1, color ="k")
savefig("Relative Displacement Spectra.png")
show()

plt.figure()
plt.title("Relative Velocity Spectra")
plt.plot(T12,lista_u1p_max ,linewidth = 1.5,label="Sv Los Gatos 3%",alpha=0.7,color="b")
plt.plot(T12,lista_u1p_max_pseudo ,linewidth = 1.5,label="Spv Los Gatos 3%",alpha=0.7,color="k",linestyle='dashed')
plt.plot(T12,listaT_u1p_max ,linewidth = 1.5,label="Sv Tarzana 3%",alpha=0.7,color="r")
plt.plot(T12,listaT_u1p_max_pseudo ,linewidth = 1.5,label="Spv Tarzana 3%",alpha=0.7,color="g",linestyle='dashed')
plt.legend(loc="upper right",fontsize="medium")
plt.grid(b=True, alpha=0.7,linestyle="dashed", color="dimgray",linewidth=0.5)
plt.xlabel("Period [s]",size = 12)
plt.ylabel("Sv , Spv [m/s]",size = 12)
plt.xlim(0,3)
plt.ylim(0,5)
# plt.axhline(0, linewidth=1, color ="k")
savefig("Relative Velocity Spectra.png")
show()


plt.figure()
plt.title("Relative Acceleration Spectra")
plt.plot(T12,lista_u2p_max ,linewidth = 1.5,label="Sa Los Gatos 3%",alpha=0.7,color="b")
plt.plot(T12,lista_u2p_max_pseudo ,linewidth = 1.5,label="Spa Los Gatos 3%",alpha=0.7,color="k",linestyle='dashed')
plt.plot(T12,listaT_u2p_max ,linewidth = 1.5,label="Sa Tarzana 3%",alpha=0.7,color="r")
plt.plot(T12,listaT_u2p_max_pseudo ,linewidth = 1.5,label="Spa Tarzana 3%",alpha=0.7,color="g",linestyle='dashed')
plt.legend(loc="upper right",fontsize="medium")
plt.grid(b=True, alpha=0.7,linestyle="dashed", color="dimgray",linewidth=0.5)
plt.xlabel("Period [s]",size = 12)
plt.ylabel("Sa , Spa [m/s2]",size = 12)
plt.xlim(0,3)
plt.ylim(0,65)
# plt.axhline(0, linewidth=1, color ="k")
savefig("Relative Acceleration Spectra.png")
show()


lista_u_max=[]
lista_u1p_max=[]
lista_u2p_max=[]

lista_u1p_max_pseudo=[]
lista_u2p_max_pseudo=[]

# LOS GATOS
for j in T12: 
    chi=0.05
    masa = 1.0 #kg
    w=(2*np.pi)/j
    c= 2*chi*w #kg/s 
    k=w*w*masa #kg/s2
    dt1=1/50
    
    # Condiciones iniciales
    
    u0=0
    u1p0=0
    u2p0= G2[0]/-masa

    lista_u=[u0]
    lista_u1p=[u1p0]
    lista_u2p=[u2p0]
    lista_p_g=[]
    
    a1=(1/(beta*(dt1*dt1))*masa)+(gamma/(beta*dt1))*c
    a2=(1/(beta*dt1))*masa + (gamma/beta-1)*c
    a3=((1/(2*beta))-1)*masa + dt1*(gamma/(2*beta)-1)*c

    k_g= k + a1

    for i in range(len(G2)-1):
     
        p_g=G2[i+1]+a1*lista_u[i]+a2*lista_u1p[i]+a3*lista_u2p[i]
        lista_p_g.append(p_g)
        
        u_i=lista_p_g[i]/k_g
        
        u1p_i=(gamma/(beta*dt1))*(u_i-lista_u[i])+(1-(gamma/beta))*lista_u1p[i]+dt1*(1-(gamma/(2*beta)))*lista_u2p[i]

        u2p_i=(1/(beta*(dt1*dt1)))*(u_i-lista_u[i])-(1/(beta*dt1))*lista_u1p[i]-(1/(2*beta)-1)*lista_u2p[i]

        lista_u.append(u_i)
        lista_u1p.append(u1p_i)
        lista_u2p.append(u2p_i)
        
    lista_u_max.append(abs(max(lista_u,key=abs)))
    lista_u1p_max.append(abs(max(lista_u1p,key=abs)))   
    lista_u2p_max.append(abs(max(lista_u2p,key=abs)))
    
    lista_u1p_max_pseudo.append(w*abs(max(lista_u,key=abs)))
    lista_u2p_max_pseudo.append(w*w*abs(max(lista_u,key=abs))) 

# TARZANA

listaT_u_max=[]
listaT_u1p_max=[]
listaT_u2p_max=[]

listaT_u1p_max_pseudo=[]
listaT_u2p_max_pseudo=[]

for j in T12: 
    chi=0.05
    masa = 1.0 #kg
    w=(2*np.pi)/j
    c= 2*chi*w #kg/s 
    k=w*w*masa #kg/s2
    dt1=1/50
    
    # Condiciones iniciales
    
    u0=0
    u1p0=0
    u2p0= T2[0]/-masa

    lista_u=[u0]
    lista_u1p=[u1p0]
    lista_u2p=[u2p0]
    lista_p_g=[]
    
    a1=(1/(beta*(dt1*dt1))*masa)+(gamma/(beta*dt1))*c
    a2=(1/(beta*dt1))*masa + (gamma/beta-1)*c
    a3=((1/(2*beta))-1)*masa + dt1*(gamma/(2*beta)-1)*c

    k_g= k + a1

    for i in range(len(T2)-1):
     
        p_g=T2[i+1]+a1*lista_u[i]+a2*lista_u1p[i]+a3*lista_u2p[i]
        lista_p_g.append(p_g)
        
        u_i=lista_p_g[i]/k_g
        
        u1p_i=(gamma/(beta*dt1))*(u_i-lista_u[i])+(1-(gamma/beta))*lista_u1p[i]+dt1*(1-(gamma/(2*beta)))*lista_u2p[i]

        u2p_i=(1/(beta*(dt1*dt1)))*(u_i-lista_u[i])-(1/(beta*dt1))*lista_u1p[i]-(1/(2*beta)-1)*lista_u2p[i]

        lista_u.append(u_i)
        lista_u1p.append(u1p_i)
        lista_u2p.append(u2p_i)
        
    listaT_u_max.append(abs(max(lista_u,key=abs)))
    listaT_u1p_max.append(abs(max(lista_u1p,key=abs)))
    listaT_u2p_max.append(abs(max(lista_u2p,key=abs)))    
    
    listaT_u1p_max_pseudo.append(w*abs(max(lista_u,key=abs)))
    listaT_u2p_max_pseudo.append(w*w*abs(max(lista_u,key=abs))) 
    
    
plt.figure()
plt.title("Relative Displacement Spectra")
plt.plot(T12,lista_u_max ,linewidth = 1.5,label="Los Gatos 5%",alpha=0.7,color="b")
plt.plot(T12,listaT_u_max ,linewidth = 1.5,label="Tarzana 5%",alpha=0.7,color="r")
plt.legend(loc="upper right",fontsize="medium")
plt.grid(b=True, alpha=0.7,linestyle="dashed", color="dimgray",linewidth=0.5)
plt.xlabel("Period [s]",size = 12)
plt.ylabel("Sd [m]",size = 12)
plt.xlim(0,3)
plt.ylim(0,0.9)
# plt.axhline(0, linewidth=1, color ="k")
savefig("Relative Displacement Spectra 5%.png")
show()

plt.figure()
plt.title("Relative Velocity Spectra")
plt.plot(T12,lista_u1p_max ,linewidth = 1.5,label="Sv Los Gatos 5%",alpha=0.7,color="b")
plt.plot(T12,lista_u1p_max_pseudo ,linewidth = 1.5,label="Spv Los Gatos 5%",alpha=0.7,color="k",linestyle='dashed')
plt.plot(T12,listaT_u1p_max ,linewidth = 1.5,label="Sv Tarzana 5%",alpha=0.7,color="r")
plt.plot(T12,listaT_u1p_max_pseudo ,linewidth = 1.5,label="Spv Tarzana 5%",alpha=0.7,color="g",linestyle='dashed')
plt.legend(loc="upper right",fontsize="medium")
plt.grid(b=True, alpha=0.7,linestyle="dashed", color="dimgray",linewidth=0.5)
plt.xlabel("Period [s]",size = 12)
plt.ylabel("Sv , Spv [m/s]",size = 12)
plt.xlim(0,3)
plt.ylim(0,5)
# plt.axhline(0, linewidth=1, color ="k")
savefig("Relative Velocity Spectra 5%.png")
show()


plt.figure()
plt.title("Relative Acceleration Spectra")
plt.plot(T12,lista_u2p_max ,linewidth = 1.5,label="Sa Los Gatos 5%",alpha=0.7,color="b")
plt.plot(T12,lista_u2p_max_pseudo ,linewidth = 1.5,label="Spa Los Gatos 5%",alpha=0.7,color="k",linestyle='dashed')
plt.plot(T12,listaT_u2p_max ,linewidth = 1.5,label="Sa Tarzana 5%",alpha=0.7,color="r")
plt.plot(T12,listaT_u2p_max_pseudo ,linewidth = 1.5,label="Spa Tarzana 5%",alpha=0.7,color="g",linestyle='dashed')
plt.legend(loc="upper right",fontsize="medium")
plt.grid(b=True, alpha=0.7,linestyle="dashed", color="dimgray",linewidth=0.5)
plt.xlabel("Period [s]",size = 12)
plt.ylabel("Sa , Spa [m/s2]",size = 12)
plt.xlim(0,3)
plt.ylim(0,65)
# plt.axhline(0, linewidth=1, color ="k")
savefig("Relative Acceleration Spectra 5%.png")
show()


# LOS GATOS PREGUNTA 3, SISTEMA 1 y 2 chi =0.0312

lista_u_max=[]
lista_u1p_max=[]
lista_u2p_max=[]

lista_u1p_max_pseudo=[]
lista_u2p_max_pseudo=[]
for j in T12: 
    chi=0.0312
    masa = 1.0 #kg
    w=(2*np.pi)/j
    c= 2*chi*w #kg/s 
    k=w*w*masa #kg/s2
    dt1=1/50
    
    # Condiciones iniciales
    
    u0=0
    u1p0=0
    u2p0= G2[0]/-masa

    lista_u=[u0]
    lista_u1p=[u1p0]
    lista_u2p=[u2p0]
    lista_p_g=[]
    
    a1=(1/(beta*(dt1*dt1))*masa)+(gamma/(beta*dt1))*c
    a2=(1/(beta*dt1))*masa + (gamma/beta-1)*c
    a3=((1/(2*beta))-1)*masa + dt1*(gamma/(2*beta)-1)*c

    k_g= k + a1

    for i in range(len(G2)-1):
     
        p_g=G2[i+1]+a1*lista_u[i]+a2*lista_u1p[i]+a3*lista_u2p[i]
        lista_p_g.append(p_g)
        
        u_i=lista_p_g[i]/k_g
        
        u1p_i=(gamma/(beta*dt1))*(u_i-lista_u[i])+(1-(gamma/beta))*lista_u1p[i]+dt1*(1-(gamma/(2*beta)))*lista_u2p[i]

        u2p_i=(1/(beta*(dt1*dt1)))*(u_i-lista_u[i])-(1/(beta*dt1))*lista_u1p[i]-(1/(2*beta)-1)*lista_u2p[i]

        lista_u.append(u_i)
        lista_u1p.append(u1p_i)
        lista_u2p.append(u2p_i)
        
    lista_u_max.append(abs(max(lista_u,key=abs)))
    lista_u1p_max.append(abs(max(lista_u1p,key=abs)))   
    lista_u2p_max.append(abs(max(lista_u2p,key=abs)))
    
    lista_u1p_max_pseudo.append(w*abs(max(lista_u,key=abs)))
    lista_u2p_max_pseudo.append(w*w*abs(max(lista_u,key=abs))) 

# TARZANA

listaT_u_max=[]
listaT_u1p_max=[]
listaT_u2p_max=[]

listaT_u1p_max_pseudo=[]
listaT_u2p_max_pseudo=[]

for j in T12: 
    chi=0.0312
    masa = 1.0 #kg
    w=(2*np.pi)/j
    c= 2*chi*w #kg/s 
    k=w*w*masa #kg/s2
    dt1=1/50
    
    # Condiciones iniciales
    
    u0=0
    u1p0=0
    u2p0= T2[0]/-masa

    lista_u=[u0]
    lista_u1p=[u1p0]
    lista_u2p=[u2p0]
    lista_p_g=[]
    
    a1=(1/(beta*(dt1*dt1))*masa)+(gamma/(beta*dt1))*c
    a2=(1/(beta*dt1))*masa + (gamma/beta-1)*c
    a3=((1/(2*beta))-1)*masa + dt1*(gamma/(2*beta)-1)*c

    k_g= k + a1

    for i in range(len(T2)-1):
     
        p_g=T2[i+1]+a1*lista_u[i]+a2*lista_u1p[i]+a3*lista_u2p[i]
        lista_p_g.append(p_g)
        
        u_i=lista_p_g[i]/k_g
        
        u1p_i=(gamma/(beta*dt1))*(u_i-lista_u[i])+(1-(gamma/beta))*lista_u1p[i]+dt1*(1-(gamma/(2*beta)))*lista_u2p[i]

        u2p_i=(1/(beta*(dt1*dt1)))*(u_i-lista_u[i])-(1/(beta*dt1))*lista_u1p[i]-(1/(2*beta)-1)*lista_u2p[i]

        lista_u.append(u_i)
        lista_u1p.append(u1p_i)
        lista_u2p.append(u2p_i)
        
    listaT_u_max.append(abs(max(lista_u,key=abs)))
    listaT_u1p_max.append(abs(max(lista_u1p,key=abs)))
    listaT_u2p_max.append(abs(max(lista_u2p,key=abs)))    
    
    listaT_u1p_max_pseudo.append(w*abs(max(lista_u,key=abs)))
    listaT_u2p_max_pseudo.append(w*w*abs(max(lista_u,key=abs))) 
    
    
plt.figure()
plt.title("Relative Displacement Spectra System 1 & 2")
plt.plot(T12,lista_u_max ,linewidth = 1.5,label="Los Gatos 3.12%",alpha=0.7,color="b")
plt.plot(T12,listaT_u_max ,linewidth = 1.5,label="Tarzana 3.12%",alpha=0.7,color="r")
plt.legend(loc="upper right",fontsize="medium")
plt.grid(b=True, alpha=0.7,linestyle="dashed", color="dimgray",linewidth=0.5)
plt.xlabel("Period [s]",size = 12)
plt.ylabel("Sd [m]",size = 12)
plt.xlim(0,3)
plt.ylim(0,0.9)
# plt.axhline(0, linewidth=1, color ="k")
savefig("PseudoDispQ3a.png")
show()

plt.figure()
plt.title("Relative Velocity Pseudo-Spectra System 1 & 2")
# plt.plot(T12,lista_u1p_max ,linewidth = 1.5,label="Sv Los Gatos 3%",alpha=0.7,color="b")
plt.plot(T12,lista_u1p_max_pseudo ,linewidth = 1.5,label="Spv Los Gatos 3.12%",alpha=0.7,color="k",linestyle='dashed')
# plt.plot(T12,listaT_u1p_max ,linewidth = 1.5,label="Sv Tarzana 3%",alpha=0.7,color="r")
plt.plot(T12,listaT_u1p_max_pseudo ,linewidth = 1.5,label="Spv Tarzana 3.12%",alpha=0.7,color="g",linestyle='dashed')
plt.legend(loc="upper right",fontsize="medium")
plt.grid(b=True, alpha=0.7,linestyle="dashed", color="dimgray",linewidth=0.5)
plt.xlabel("Period [s]",size = 12)
plt.ylabel("Sv , Spv [m/s]",size = 12)
plt.xlim(0,3)
plt.ylim(0,5)
# plt.axhline(0, linewidth=1, color ="k")
savefig("PseudoVelQ3a.png")
show()


plt.figure()
plt.title("Relative Acceleration Pseudo-Spectra System 1 & 2")
# plt.plot(T12,lista_u2p_max ,linewidth = 1.5,label="Sa Los Gatos 3.12%",alpha=0.7,color="b")
plt.plot(T12,lista_u2p_max_pseudo ,linewidth = 1.5,label="Spa Los Gatos 3.12%",alpha=0.7,color="k",linestyle='dashed')
# plt.plot(T12,listaT_u2p_max ,linewidth = 1.5,label="Sa Tarzana 3%",alpha=0.7,color="r")
plt.plot(T12,listaT_u2p_max_pseudo ,linewidth = 1.5,label="Spa Tarzana 3.12%",alpha=0.7,color="g",linestyle='dashed')
plt.legend(loc="upper right",fontsize="medium")
plt.grid(b=True, alpha=0.7,linestyle="dashed", color="dimgray",linewidth=0.5)
plt.xlabel("Period [s]",size = 12)
plt.ylabel("Sa , Spa [m/s2]",size = 12)
plt.xlim(0,3)
plt.ylim(0,65)
# plt.axhline(0, linewidth=1, color ="k")
savefig("PseudoAcceQ3a.png")
show()


# LOS GATOS PREGUNTA 3, SISTEMA 1 y 2 chi =0.05

lista_u_max=[]
lista_u1p_max=[]
lista_u2p_max=[]

lista_u1p_max_pseudo=[]
lista_u2p_max_pseudo=[]
for j in T12: 
    chi=0.05
    masa = 1.0 #kg
    w=(2*np.pi)/j
    c= 2*chi*w #kg/s 
    k=w*w*masa #kg/s2
    dt1=1/50
    
    # Condiciones iniciales
    
    u0=0
    u1p0=0
    u2p0= G2[0]/-masa

    lista_u=[u0]
    lista_u1p=[u1p0]
    lista_u2p=[u2p0]
    lista_p_g=[]
    
    a1=(1/(beta*(dt1*dt1))*masa)+(gamma/(beta*dt1))*c
    a2=(1/(beta*dt1))*masa + (gamma/beta-1)*c
    a3=((1/(2*beta))-1)*masa + dt1*(gamma/(2*beta)-1)*c

    k_g= k + a1

    for i in range(len(G2)-1):
     
        p_g=G2[i+1]+a1*lista_u[i]+a2*lista_u1p[i]+a3*lista_u2p[i]
        lista_p_g.append(p_g)
        
        u_i=lista_p_g[i]/k_g
        
        u1p_i=(gamma/(beta*dt1))*(u_i-lista_u[i])+(1-(gamma/beta))*lista_u1p[i]+dt1*(1-(gamma/(2*beta)))*lista_u2p[i]

        u2p_i=(1/(beta*(dt1*dt1)))*(u_i-lista_u[i])-(1/(beta*dt1))*lista_u1p[i]-(1/(2*beta)-1)*lista_u2p[i]

        lista_u.append(u_i)
        lista_u1p.append(u1p_i)
        lista_u2p.append(u2p_i)
        
    lista_u_max.append(abs(max(lista_u,key=abs)))
    lista_u1p_max.append(abs(max(lista_u1p,key=abs)))   
    lista_u2p_max.append(abs(max(lista_u2p,key=abs)))
    
    lista_u1p_max_pseudo.append(w*abs(max(lista_u,key=abs)))
    lista_u2p_max_pseudo.append(w*w*abs(max(lista_u,key=abs))) 

# TARZANA

listaT_u_max=[]
listaT_u1p_max=[]
listaT_u2p_max=[]

listaT_u1p_max_pseudo=[]
listaT_u2p_max_pseudo=[]

for j in T12: 
    chi=0.05
    masa = 1.0 #kg
    w=(2*np.pi)/j
    c= 2*chi*w #kg/s 
    k=w*w*masa #kg/s2
    dt1=1/50
    
    # Condiciones iniciales
    
    u0=0
    u1p0=0
    u2p0= T2[0]/-masa

    lista_u=[u0]
    lista_u1p=[u1p0]
    lista_u2p=[u2p0]
    lista_p_g=[]
    
    a1=(1/(beta*(dt1*dt1))*masa)+(gamma/(beta*dt1))*c
    a2=(1/(beta*dt1))*masa + (gamma/beta-1)*c
    a3=((1/(2*beta))-1)*masa + dt1*(gamma/(2*beta)-1)*c

    k_g= k + a1

    for i in range(len(T2)-1):
     
        p_g=T2[i+1]+a1*lista_u[i]+a2*lista_u1p[i]+a3*lista_u2p[i]
        lista_p_g.append(p_g)
        
        u_i=lista_p_g[i]/k_g
        
        u1p_i=(gamma/(beta*dt1))*(u_i-lista_u[i])+(1-(gamma/beta))*lista_u1p[i]+dt1*(1-(gamma/(2*beta)))*lista_u2p[i]

        u2p_i=(1/(beta*(dt1*dt1)))*(u_i-lista_u[i])-(1/(beta*dt1))*lista_u1p[i]-(1/(2*beta)-1)*lista_u2p[i]

        lista_u.append(u_i)
        lista_u1p.append(u1p_i)
        lista_u2p.append(u2p_i)
        
    listaT_u_max.append(abs(max(lista_u,key=abs)))
    listaT_u1p_max.append(abs(max(lista_u1p,key=abs)))
    listaT_u2p_max.append(abs(max(lista_u2p,key=abs)))    
    
    listaT_u1p_max_pseudo.append(w*abs(max(lista_u,key=abs)))
    listaT_u2p_max_pseudo.append(w*w*abs(max(lista_u,key=abs))) 
    
    
plt.figure()
plt.title("Relative Displacement Spectra")
plt.plot(T12,lista_u_max ,linewidth = 1.5,label="Los Gatos 5%",alpha=0.7,color="b")
plt.plot(T12,listaT_u_max ,linewidth = 1.5,label="Tarzana 5%",alpha=0.7,color="r")
plt.legend(loc="upper right",fontsize="medium")
plt.grid(b=True, alpha=0.7,linestyle="dashed", color="dimgray",linewidth=0.5)
plt.xlabel("Period [s]",size = 12)
plt.ylabel("Sd [m]",size = 12)
plt.xlim(0,3)
plt.ylim(0,0.9)
# plt.axhline(0, linewidth=1, color ="k")
savefig("PseudoDispQ3b.png")
show()

plt.figure()
plt.title("Relative Velocity Pseudo-Spectra")
# plt.plot(T12,lista_u1p_max ,linewidth = 1.5,label="Sv Los Gatos 3%",alpha=0.7,color="b")
plt.plot(T12,lista_u1p_max_pseudo ,linewidth = 1.5,label="Spv Los Gatos 5%",alpha=0.7,color="k",linestyle='dashed')
# plt.plot(T12,listaT_u1p_max ,linewidth = 1.5,label="Sv Tarzana 3%",alpha=0.7,color="r")
plt.plot(T12,listaT_u1p_max_pseudo ,linewidth = 1.5,label="Spv Tarzana 5%",alpha=0.7,color="g",linestyle='dashed')
plt.legend(loc="upper right",fontsize="medium")
plt.grid(b=True, alpha=0.7,linestyle="dashed", color="dimgray",linewidth=0.5)
plt.xlabel("Period [s]",size = 12)
plt.ylabel("Sv , Spv [m/s]",size = 12)
plt.xlim(0,3)
plt.ylim(0,5)
# plt.axhline(0, linewidth=1, color ="k")
savefig("PseudoVelQ3b.png")
show()


plt.figure()
plt.title("Relative Acceleration Pseudo-Spectra")
# plt.plot(T12,lista_u2p_max ,linewidth = 1.5,label="Sa Los Gatos 3.12%",alpha=0.7,color="b")
plt.plot(T12,lista_u2p_max_pseudo ,linewidth = 1.5,label="Spa Los Gatos 5%",alpha=0.7,color="k",linestyle='dashed')
# plt.plot(T12,listaT_u2p_max ,linewidth = 1.5,label="Sa Tarzana 3%",alpha=0.7,color="r")
plt.plot(T12,listaT_u2p_max_pseudo ,linewidth = 1.5,label="Spa Tarzana 5%",alpha=0.7,color="g",linestyle='dashed')
plt.legend(loc="upper right",fontsize="medium")
plt.grid(b=True, alpha=0.7,linestyle="dashed", color="dimgray",linewidth=0.5)
plt.xlabel("Period [s]",size = 12)
plt.ylabel("Sa , Spa [m/s2]",size = 12)
plt.xlim(0,3)
plt.ylim(0,65)
# plt.axhline(0, linewidth=1, color ="k")
savefig("PseudoAcceQ3b.png")
show()

# #//////////////////////////////////////////////////////////////////////////////

# # #/// [2] //////////////////////////////////////////////////////////////////////


# # # ### SISTEMA 1 ###############################################################

# # # --- LOS GATOS ---------------------------------------------------------------
# # #Newmark’s method: linear acceleration.

masa = 1500. #kg
beta=1/6
gamma=1/2
delta_t=dt #s
c= 490.333 #kg/s
k=41187.93 #kg/s2


"Condiciones Iniciales"

u0=0

u1p0=0

u2p0= G2[0]/-masa


a1=(1/(beta*(delta_t*delta_t))*masa)+(gamma/(beta*delta_t))*c

a2=(1/(beta*delta_t))*masa + (gamma/beta-1)*c

a3=((1/(2*beta))-1)*masa + delta_t*(gamma/(2*beta)-1)*c

k_g= k + a1

"2.0""Calculos para pasos de tiempo i=0,1,2,....n"


lista_u=[u0]

lista_u1p=[u1p0]

lista_u2p=[u2p0]

lista_p_g=[]
ct=0
for i in G2:
    if ct+1 == len(G2):
        break
    p_g=-masa*G2[ct+1]+a1*lista_u[ct]+a2*lista_u1p[ct]+a3*lista_u2p[ct]
    lista_p_g.append(p_g)
    
    u_i=lista_p_g[ct]/k_g
    
    u1p_i=(gamma/(beta*delta_t))*(u_i-lista_u[ct])+(1-(gamma/beta))*lista_u1p[ct]+delta_t*(1-(gamma/(2*beta)))*lista_u2p[ct]

    u2p_i=(1/(beta*(delta_t*delta_t)))*(u_i-lista_u[ct])-(1/(beta*delta_t))*lista_u1p[ct]-(1/(2*beta)-1)*lista_u2p[ct]
    
    lista_u.append(u_i)
    lista_u1p.append(u1p_i)
    lista_u2p.append(u2p_i)
    ct+=1
     
#Central difference method
       
k_gorro =masa/(dt*dt)+c/(2*dt)

u_cd = [0]*2000        
udot = [0]*2000      
u2dot = [0]*2000     
rgor = [0]*2000   

u_cd[0] = 0
udot[0] = 0
u2dot[0] = G2[0]/masa

for i in range(0,1999):
    if i == 0:
        u_cd[i+1] = dt*dt*u2dot[i]/2
        udot[i]=0
        u2dot[i]= (u_cd[i+1]-2*u_cd[i]+u2dot[i]*((dt*dt)/2))/(dt*dt)
        
    else:
        rgor[i]=-masa*G2[i]-((masa/(dt*dt)-c/(2*dt))*u_cd[i-1]) - (k - (2*masa)/(dt*dt))*u_cd[i]
        u_cd[i+1]=rgor[i]/k_gorro
        udot[i]= (u_cd[i+1]-u_cd[i-1])/(2*dt)
        u2dot[i]=(u_cd[i+1]-2*u_cd[i]+u_cd[i-1])/(dt*dt)
       
       
lista_u_cd= u_cd 
lista_u1p_cd= udot
lista_u2p_cd= u2dot

#DUHAMELS INTEGRAL AND RECTANGULAR SUM

w = np.sqrt(k/masa)
Periodo = 2*np.pi/w
xi = c/(2*w*masa)

wd=w*np.sqrt(1-xi*xi)
N = len(G2)

funA=[0]*2000
funB=[0]*2000
u_di=[0]*2000
u_pdi=[0]*2000
u_2pi=[0]*2000

A=[0]*2000
B=[0]*2000
for i in range(0,N):
    tau=dt*(i-1)
    funA[i]= np.exp(xi*w*tau)*-masa*G2[i]*np.cos(wd*tau)
    funB[i]= np.exp(xi*w*tau)*-masa*G2[i]*np.sin(wd*tau)

A[0]=0
B[0]=0
u_di[0]=0
u_pdi[0]=0
u_2pi[0]=0

for i in range(1,N-1):
    A[i]=A[i-1]+1/(masa*wd)*dt*funA[i]
    B[i]=B[i-1]+1/(masa*wd)*dt*funB[i]
    u_di[i]=A[i]*np.exp(-xi*w*dt*i)*np.sin(wd*dt*i)-B[i]*np.exp(-xi*w*dt*i)*np.cos(wd*dt*i)
    #
    u_pdi[i]=(-(xi*A[i]*w)+wd*B[i])*np.exp(-xi*dt*i*w)*np.sin(wd*dt*i)+(xi*B[i]*w+wd*A[i])*np.cos(wd*dt*i)*np.exp(-xi*dt*i*w)
    u_2pi[i]=(xi*xi*A[i]*w*w-2*wd*xi*B[i]*w-wd*wd*A[i])*np.exp(-xi*dt*i*w)*np.sin(wd*dt*i)+(wd*wd*B[i]-(xi*xi*B[i]*w*w+2*wd*xi*A[i]*w))*np.cos(wd*dt*i)*np.exp(-xi*dt*i*w)-G2[i]
# -----------------------------------------------------------------------------

#
### PLOTEO LOS GATOS ###   
#
 
plt.figure()
plt.title("Displacement System 1 Los Gatos")
plt.plot(tg,lista_u ,linewidth = 1.5,label="Newmark",alpha=0.7,color="b")
plt.plot(tg,u_di ,linewidth = 1.5,label="Duhamel's",alpha=0.7,color="g",linestyle="dashed")
plt.plot(tg,lista_u_cd ,linestyle="dotted",linewidth = 1.5,label="Central difference",alpha=0.7,color="r")
plt.legend(loc="upper right",fontsize="larger")
plt.grid(b=True, alpha=0.7,linestyle="dashed", color="dimgray",linewidth=0.5)
plt.xlabel("Time [s]",size = 12)
plt.ylabel("Displacement [m]",size = 12)
plt.xlim(0,40)
plt.axhline(0, linewidth=1, color ="k")
savefig("Earthquake_Record_displ_LG_S1.png")
show()

plt.figure()
plt.title("Velocity System 1 Los Gatos")
plt.plot(tg,lista_u1p ,linewidth = 1.5,label="Newmark",alpha=0.7,color = "k")
plt.plot(tg,u_pdi ,linewidth = 1.5,label="Duhamel's",alpha=0.7,color="g",linestyle="dashed")
plt.plot(tg,lista_u1p_cd,linestyle="dotted",linewidth = 1.5,label="Central difference",alpha=0.7,color="y")
plt.legend(loc="upper right",fontsize="larger")
plt.grid(b=True, alpha=0.7,linestyle="dashed", color="dimgray",linewidth=0.5)
plt.xlabel("Time [s]",size = 12)
plt.ylabel("velocity [m/s]",size = 12)
plt.xlim(0,40)
plt.axhline(0, linewidth=1, color ="k")
savefig("Earthquake_Record_vel_LG_S1.png")
show()

plt.figure()
plt.title("Acceleration System 1 Los Gatos")
plt.plot(tg,lista_u2p ,linewidth = 1.5,label="Newmark",alpha=0.7)
plt.plot(tg,u_2pi ,linewidth = 1.5,label="Duhamel's",alpha=0.7,color="g",linestyle="dashed")
plt.plot(tg,lista_u2p_cd,linestyle="dashed",linewidth = 1.5,label="Central difference",alpha=0.7)
plt.legend(loc="upper right",fontsize="larger")
plt.grid(b=True, alpha=0.7,linestyle="dashed", color="dimgray",linewidth=0.5)
plt.xlabel("Time [s]",size = 12)
plt.ylabel("Acceleration [m/s2]",size = 12)
plt.xlim(0,40)
plt.axhline(0, linewidth=1, color ="k")
savefig("Earthquake_Record_Acce_LG_S1.png")
show()

# --- TARZANA -----------------------------------------------------------------

#Newmark’s method: linear acceleration.

masa = 1500. #kg
beta=1/6
gamma=1/2
delta_t=dt #s
c= 490.333 #kg/s
k=41187.93 #kg/s2


"Condiciones Iniciales"

u0=0

u1p0=0

u2p0= T2[0]/-masa


a1=(1/(beta*(delta_t*delta_t))*masa)+(gamma/(beta*delta_t))*c

a2=(1/(beta*delta_t))*masa + (gamma/beta-1)*c

a3=((1/(2*beta))-1)*masa + delta_t*(gamma/(2*beta)-1)*c

k_g= k + a1

"2.0""Calculos para pasos de tiempo i=0,1,2,....n"


lista_u=[u0]

lista_u1p=[u1p0]

lista_u2p=[u2p0]

lista_p_g=[]
ct=0
for i in T2:
    if ct+1 == len(T2):
        break
    p_g=-masa*T2[ct+1]+a1*lista_u[ct]+a2*lista_u1p[ct]+a3*lista_u2p[ct]
    lista_p_g.append(p_g)
    
    u_i=lista_p_g[ct]/k_g
    
    u1p_i=(gamma/(beta*delta_t))*(u_i-lista_u[ct])+(1-(gamma/beta))*lista_u1p[ct]+delta_t*(1-(gamma/(2*beta)))*lista_u2p[ct]

    u2p_i=(1/(beta*(delta_t*delta_t)))*(u_i-lista_u[ct])-(1/(beta*delta_t))*lista_u1p[ct]-(1/(2*beta)-1)*lista_u2p[ct]
    
    lista_u.append(u_i)
    lista_u1p.append(u1p_i)
    lista_u2p.append(u2p_i)
    ct+=1
     
#Central difference method
       
k_gorro =masa/(dt*dt)+c/(2*dt)

u_cd = [0]*len(T2)       
udot = [0]*len(T2)     
u2dot = [0]*len(T2)     
rgor = [0]*len(T2)  

u_cd[0] = 0
udot[0] = 0
u2dot[0] = T2[0]/masa

for i in range(0,len(T2)-1):
    if i == 0:
        u_cd[i+1] = dt*dt*u2dot[i]/2
        udot[i]=0
        u2dot[i]= (u_cd[i+1]-2*u_cd[i]+u2dot[i]*((dt*dt)/2))/(dt*dt)
        
    else:
        rgor[i]=-masa*T2[i]-((masa/(dt*dt)-c/(2*dt))*u_cd[i-1]) - (k - (2*masa)/(dt*dt))*u_cd[i]
        u_cd[i+1]=rgor[i]/k_gorro
        udot[i]= (u_cd[i+1]-u_cd[i-1])/(2*dt)
        u2dot[i]=(u_cd[i+1]-2*u_cd[i]+u_cd[i-1])/(dt*dt)
       
       
lista_u_cd= u_cd 
lista_u1p_cd= udot
lista_u2p_cd= u2dot

#DUHAMELS INTEGRAL AND RECTANGULAR SUM

w = np.sqrt(k/masa)
Periodo = 2*np.pi/w
xi = c/(2*w*masa)

wd=w*np.sqrt(1-xi*xi)
N = len(T2)

funA=[0]*3001
funB=[0]*3001
u_di=[0]*3001
u_pdi=[0]*3001
u_2pi=[0]*3001

A=[0]*3001
B=[0]*3001
for i in range(0,N):
    tau=dt*(i-1)
    funA[i]= np.exp(xi*w*tau)*-masa*T2[i]*np.cos(wd*tau)
    funB[i]= np.exp(xi*w*tau)*-masa*T2[i]*np.sin(wd*tau)

A[0]=0
B[0]=0
u_di[0]=0
u_pdi[0]=0
u_2pi[0]=0

for i in range(1,N-1):
    A[i]=A[i-1]+1/(masa*wd)*dt*funA[i]
    B[i]=B[i-1]+1/(masa*wd)*dt*funB[i]
    u_di[i]=A[i]*np.exp(-xi*w*dt*i)*np.sin(wd*dt*i)-B[i]*np.exp(-xi*w*dt*i)*np.cos(wd*dt*i)
    #
    u_pdi[i]=(-(xi*A[i]*w)+wd*B[i])*np.exp(-xi*dt*i*w)*np.sin(wd*dt*i)+(xi*B[i]*w+wd*A[i])*np.cos(wd*dt*i)*np.exp(-xi*dt*i*w)
    u_2pi[i]=(xi*xi*A[i]*w*w-2*wd*xi*B[i]*w-wd*wd*A[i])*np.exp(-xi*dt*i*w)*np.sin(wd*dt*i)+(wd*wd*B[i]-(xi*xi*B[i]*w*w+2*wd*xi*A[i]*w))*np.cos(wd*dt*i)*np.exp(-xi*dt*i*w)-T2[i]
# -----------------------------------------------------------------------------

#
### PLOTEO TARZANA ###  
#

plt.figure()
plt.title("Displacement System 1 Tarzana")
plt.plot(tt,lista_u ,linewidth = 1.5,label="Newmark",alpha=0.7,color="b")
plt.plot(tt,u_di ,linewidth = 1.5,label="Duhamel's",alpha=0.7,color="g",linestyle="dashed")
plt.plot(tt,lista_u_cd ,linestyle="dotted",linewidth = 1.5,label="Central difference",alpha=0.7,color="r")
plt.legend(loc="upper right",fontsize="larger")
plt.grid(b=True, alpha=0.7,linestyle="dashed", color="dimgray",linewidth=0.5)
plt.xlabel("Time [s]",size = 12)
plt.ylabel("Displacement [m]",size = 12)
plt.xlim(0,40)
plt.axhline(0, linewidth=1, color ="k")
savefig("Earthquake_Record_displ_TZ_S1.png")
show()

plt.figure()
plt.title("Velocity System 1 Tarzana")
plt.plot(tt,lista_u1p ,linewidth = 1.5,label="Newmark",alpha=0.7,color = "k")
plt.plot(tt,u_pdi ,linewidth = 1.5,label="Duhamel's",alpha=0.7,color="g",linestyle="dashed")
plt.plot(tt,lista_u1p_cd,linestyle="dotted",linewidth = 1.5,label="Central difference",alpha=0.7,color="y")
plt.legend(loc="upper right",fontsize="larger")
plt.grid(b=True, alpha=0.7,linestyle="dashed", color="dimgray",linewidth=0.5)
plt.xlabel("Time [s]",size = 12)
plt.ylabel("velocity [m/s]",size = 12)
plt.xlim(0,40)
plt.axhline(0, linewidth=1, color ="k")
savefig("Earthquake_Record_vel_TZ_S1.png")
show()

plt.figure()
plt.title("Acceleration System 1 Tarzana")
plt.plot(tt,lista_u2p ,linewidth = 1.5,label="Newmark",alpha=0.7)
plt.plot(tt,u_2pi ,linewidth = 1.5,label="Duhamel's",alpha=0.7,color="g",linestyle="dashed")
plt.plot(tt,lista_u2p_cd,linestyle="dashed",linewidth = 1.5,label="Central difference",alpha=0.7)
plt.legend(loc="upper right",fontsize="larger")
plt.grid(b=True, alpha=0.7,linestyle="dashed", color="dimgray",linewidth=0.5)
plt.xlabel("Time [s]",size = 12)
plt.ylabel("Acceleration [m/s2]",size = 12)
plt.xlim(0,40)
plt.axhline(0, linewidth=1, color ="k")
savefig("Earthquake_Record_Acce_TZ_S1.png")
show()
# #############################################################################


# ### SISTEMA 2 ###############################################################

# --- LOS GATOS ---------------------------------------------------------------
#Newmark’s method: linear acceleration.

masa = 130. #kg
beta=1/6
gamma=1/2
delta_t=dt #s
c= 147.1 #kg/s
k=41187.93 #kg/s2


"Condiciones Iniciales"

u0=0

u1p0=0

u2p0= G2[0]/-masa


a1=(1/(beta*(delta_t*delta_t))*masa)+(gamma/(beta*delta_t))*c

a2=(1/(beta*delta_t))*masa + (gamma/beta-1)*c

a3=((1/(2*beta))-1)*masa + delta_t*(gamma/(2*beta)-1)*c

k_g= k + a1

"2.0""Calculos para pasos de tiempo i=0,1,2,....n"


lista_u=[u0]

lista_u1p=[u1p0]

lista_u2p=[u2p0]

lista_p_g=[]
ct=0
for i in G2:
    if ct+1 == len(G2):
        break
    p_g=-masa*G2[ct+1]+a1*lista_u[ct]+a2*lista_u1p[ct]+a3*lista_u2p[ct]
    lista_p_g.append(p_g)
    
    u_i=lista_p_g[ct]/k_g
    
    u1p_i=(gamma/(beta*delta_t))*(u_i-lista_u[ct])+(1-(gamma/beta))*lista_u1p[ct]+delta_t*(1-(gamma/(2*beta)))*lista_u2p[ct]

    u2p_i=(1/(beta*(delta_t*delta_t)))*(u_i-lista_u[ct])-(1/(beta*delta_t))*lista_u1p[ct]-(1/(2*beta)-1)*lista_u2p[ct]
    
    lista_u.append(u_i)
    lista_u1p.append(u1p_i)
    lista_u2p.append(u2p_i)
    ct+=1
     
#Central difference method 
       
k_gorro =masa/(dt*dt)+c/(2*dt)

u_cd = [0]*2000        
udot = [0]*2000      
u2dot = [0]*2000     
rgor = [0]*2000   

u_cd[0] = 0
udot[0] = 0
u2dot[0] = G2[0]/masa

for i in range(0,1999):
    if i == 0:
        u_cd[i+1] = dt*dt*u2dot[i]/2
        udot[i]=0
        u2dot[i]= (u_cd[i+1]-2*u_cd[i]+u2dot[i]*((dt*dt)/2))/(dt*dt)
        
    else:
        rgor[i]=-masa*G2[i]-((masa/(dt*dt)-c/(2*dt))*u_cd[i-1]) - (k - (2*masa)/(dt*dt))*u_cd[i]
        u_cd[i+1]=rgor[i]/k_gorro
        udot[i]= (u_cd[i+1]-u_cd[i-1])/(2*dt)
        u2dot[i]=(u_cd[i+1]-2*u_cd[i]+u_cd[i-1])/(dt*dt)
       
       
lista_u_cd= u_cd 
lista_u1p_cd= udot
lista_u2p_cd= u2dot

#DUHAMELS INTEGRAL AND RECTANGULAR SUM

w = np.sqrt(k/masa)
Periodo = 2*np.pi/w
xi = c/(2*w*masa)

wd=w*np.sqrt(1-xi*xi)
N = len(G2)

funA=[0]*2000
funB=[0]*2000
u_di=[0]*2000
u_pdi=[0]*2000
u_2pi=[0]*2000

A=[0]*2000
B=[0]*2000
for i in range(0,N):
    tau=dt*(i-1)
    funA[i]= np.exp(xi*w*tau)*-masa*G2[i]*np.cos(wd*tau)
    funB[i]= np.exp(xi*w*tau)*-masa*G2[i]*np.sin(wd*tau)

A[0]=0
B[0]=0
u_di[0]=0
u_pdi[0]=0
u_2pi[0]=0

for i in range(1,N-1):
    A[i]=A[i-1]+1/(masa*wd)*dt*funA[i]
    B[i]=B[i-1]+1/(masa*wd)*dt*funB[i]
    u_di[i]=A[i]*np.exp(-xi*w*dt*i)*np.sin(wd*dt*i)-B[i]*np.exp(-xi*w*dt*i)*np.cos(wd*dt*i)
    #
    u_pdi[i]=(-(xi*A[i]*w)+wd*B[i])*np.exp(-xi*dt*i*w)*np.sin(wd*dt*i)+(xi*B[i]*w+wd*A[i])*np.cos(wd*dt*i)*np.exp(-xi*dt*i*w)
    u_2pi[i]=(xi*xi*A[i]*w*w-2*wd*xi*B[i]*w-wd*wd*A[i])*np.exp(-xi*dt*i*w)*np.sin(wd*dt*i)+(wd*wd*B[i]-(xi*xi*B[i]*w*w+2*wd*xi*A[i]*w))*np.cos(wd*dt*i)*np.exp(-xi*dt*i*w)-G2[i]
# -----------------------------------------------------------------------------

#
### PLOTEO LOS GATOS ###    
#

plt.figure()
plt.title("Displacement System 2 Los Gatos")
plt.plot(tg,lista_u ,linewidth = 1.5,label="Newmark",alpha=0.7,color="b")
plt.plot(tg,u_di ,linewidth = 1.5,label="Duhamel's",alpha=0.7,color="g",linestyle="dashed")
plt.plot(tg,lista_u_cd ,linestyle="dotted",linewidth = 1.5,label="Central difference",alpha=0.7,color="r")
plt.legend(loc="upper right",fontsize="larger")
plt.grid(b=True, alpha=0.7,linestyle="dashed", color="dimgray",linewidth=0.5)
plt.xlabel("Time [s]",size = 12)
plt.ylabel("Displacement [m]",size = 12)
plt.xlim(0,40)
plt.axhline(0, linewidth=1, color ="k")
savefig("Earthquake_Record_displ_LG_S2.png")
show()

plt.figure()
plt.title("Velocity System 2 Los Gatos")
plt.plot(tg,lista_u1p ,linewidth = 1.5,label="Newmark",alpha=0.7,color = "k")
plt.plot(tg,u_pdi ,linewidth = 1.5,label="Duhamel's",alpha=0.7,color="g",linestyle="dashed")
plt.plot(tg,lista_u1p_cd,linestyle="dotted",linewidth = 1.5,label="Central difference",alpha=0.7,color="y")
plt.legend(loc="upper right",fontsize="larger")
plt.grid(b=True, alpha=0.7,linestyle="dashed", color="dimgray",linewidth=0.5)
plt.xlabel("Time [s]",size = 12)
plt.ylabel("velocity [m/s]",size = 12)
plt.xlim(0,40)
plt.axhline(0, linewidth=1, color ="k")
savefig("Earthquake_Record_vel_LG_S2.png")
show()

plt.figure()
plt.title("Acceleration System 2 Los Gatos")
plt.plot(tg,lista_u2p ,linewidth = 1.5,label="Newmark",alpha=0.7)
plt.plot(tg,u_2pi ,linewidth = 1.5,label="Duhamel's",alpha=0.7,color="g",linestyle="dashed")
plt.plot(tg,lista_u2p_cd,linestyle="dashed",linewidth = 1.5,label="Central difference",alpha=0.7)
plt.legend(loc="upper right",fontsize="larger")
plt.grid(b=True, alpha=0.7,linestyle="dashed", color="dimgray",linewidth=0.5)
plt.xlabel("Time [s]",size = 12)
plt.ylabel("Acceleration [m/s2]",size = 12)
plt.xlim(0,40)
plt.axhline(0, linewidth=1, color ="k")
savefig("Earthquake_Record_Acce_LG_S2.png")
show()

# --- TARZANA -----------------------------------------------------------------
#Newmark’s method: linear acceleration.


masa = 130. #kg
beta=1/6
gamma=1/2
delta_t=dt #s
c= 147.1 #kg/s
k=41187.93 #kg/s2


"Condiciones Iniciales"

u0=0

u1p0=0

u2p0= T2[0]/-masa


a1=(1/(beta*(delta_t*delta_t))*masa)+(gamma/(beta*delta_t))*c

a2=(1/(beta*delta_t))*masa + (gamma/beta-1)*c

a3=((1/(2*beta))-1)*masa + delta_t*(gamma/(2*beta)-1)*c

k_g= k + a1

"2.0""Calculos para pasos de tiempo i=0,1,2,....n"


lista_u=[u0]

lista_u1p=[u1p0]

lista_u2p=[u2p0]

lista_p_g=[]
ct=0
for i in T2:
    if ct+1 == len(T2):
        break
    p_g=-masa*T2[ct+1]+a1*lista_u[ct]+a2*lista_u1p[ct]+a3*lista_u2p[ct]
    lista_p_g.append(p_g)
    
    u_i=lista_p_g[ct]/k_g
    
    u1p_i=(gamma/(beta*delta_t))*(u_i-lista_u[ct])+(1-(gamma/beta))*lista_u1p[ct]+delta_t*(1-(gamma/(2*beta)))*lista_u2p[ct]

    u2p_i=(1/(beta*(delta_t*delta_t)))*(u_i-lista_u[ct])-(1/(beta*delta_t))*lista_u1p[ct]-(1/(2*beta)-1)*lista_u2p[ct]
    
    lista_u.append(u_i)
    lista_u1p.append(u1p_i)
    lista_u2p.append(u2p_i)
    ct+=1
     
#Central difference method
       
k_gorro =masa/(dt*dt)+c/(2*dt)

u_cd = [0]*len(T2)       
udot = [0]*len(T2)     
u2dot = [0]*len(T2)     
rgor = [0]*len(T2)  

u_cd[0] = 0
udot[0] = 0
u2dot[0] = T2[0]/masa

for i in range(0,len(T2)-1):
    if i == 0:
        u_cd[i+1] = dt*dt*u2dot[i]/2
        udot[i]=0
        u2dot[i]= (u_cd[i+1]-2*u_cd[i]+u2dot[i]*((dt*dt)/2))/(dt*dt)
        
    else:
        rgor[i]=-masa*T2[i]-((masa/(dt*dt)-c/(2*dt))*u_cd[i-1]) - (k - (2*masa)/(dt*dt))*u_cd[i]
        u_cd[i+1]=rgor[i]/k_gorro
        udot[i]= (u_cd[i+1]-u_cd[i-1])/(2*dt)
        u2dot[i]=(u_cd[i+1]-2*u_cd[i]+u_cd[i-1])/(dt*dt)
       
       
lista_u_cd= u_cd 
lista_u1p_cd= udot
lista_u2p_cd= u2dot

#DUHAMELS INTEGRAL AND RECTANGULAR SUM

w = np.sqrt(k/masa)
Periodo = 2*np.pi/w
xi = c/(2*w*masa)

wd=w*np.sqrt(1-xi*xi)
N = len(T2)

funA=[0]*3001
funB=[0]*3001
u_di=[0]*3001
u_pdi=[0]*3001
u_2pi=[0]*3001

A=[0]*3001
B=[0]*3001
for i in range(0,N):
    tau=dt*(i-1)
    funA[i]= np.exp(xi*w*tau)*-masa*T2[i]*np.cos(wd*tau)
    funB[i]= np.exp(xi*w*tau)*-masa*T2[i]*np.sin(wd*tau)

A[0]=0
B[0]=0
u_di[0]=0
u_pdi[0]=0
u_2pi[0]=0

for i in range(1,N-1):
    A[i]=A[i-1]+1/(masa*wd)*dt*funA[i]
    B[i]=B[i-1]+1/(masa*wd)*dt*funB[i]
    u_di[i]=A[i]*np.exp(-xi*w*dt*i)*np.sin(wd*dt*i)-B[i]*np.exp(-xi*w*dt*i)*np.cos(wd*dt*i)
    #
    u_pdi[i]=(-(xi*A[i]*w)+wd*B[i])*np.exp(-xi*dt*i*w)*np.sin(wd*dt*i)+(xi*B[i]*w+wd*A[i])*np.cos(wd*dt*i)*np.exp(-xi*dt*i*w)
    u_2pi[i]=(xi*xi*A[i]*w*w-2*wd*xi*B[i]*w-wd*wd*A[i])*np.exp(-xi*dt*i*w)*np.sin(wd*dt*i)+(wd*wd*B[i]-(xi*xi*B[i]*w*w+2*wd*xi*A[i]*w))*np.cos(wd*dt*i)*np.exp(-xi*dt*i*w)-T2[i]
# -----------------------------------------------------------------------------
   
#
### PLOTEO TARZANA ###
#

plt.figure()
plt.title("Displacement System 2 Tarzana")
plt.plot(tt,lista_u ,linewidth = 1.5,label="Newmark",alpha=0.7,color="b")
plt.plot(tt,u_di ,linewidth = 1.5,label="Duhamel's",alpha=0.7,color="g",linestyle="dashed")
plt.plot(tt,lista_u_cd ,linestyle="dotted",linewidth = 1.5,label="Central difference",alpha=0.7,color="r")
plt.legend(loc="upper right",fontsize="larger")
plt.grid(b=True, alpha=0.7,linestyle="dashed", color="dimgray",linewidth=0.5)
plt.xlabel("Time [s]",size = 12)
plt.ylabel("Displacement [m]",size = 12)
plt.xlim(0,40)
plt.axhline(0, linewidth=1, color ="k")
savefig("Earthquake_Record_displ_TZ_S2.png")
show()

plt.figure()
plt.title("Velocity System 2 Tarzana")
plt.plot(tt,lista_u1p ,linewidth = 1.5,label="Newmark",alpha=0.7,color = "k")
plt.plot(tt,u_pdi ,linewidth = 1.5,label="Duhamel's",alpha=0.7,color="g",linestyle="dashed")
plt.plot(tt,lista_u1p_cd,linestyle="dotted",linewidth = 1.5,label="Central difference",alpha=0.7,color="y")
plt.legend(loc="upper right",fontsize="larger")
plt.grid(b=True, alpha=0.7,linestyle="dashed", color="dimgray",linewidth=0.5)
plt.xlabel("Time [s]",size = 12)
plt.ylabel("velocity [m/s]",size = 12)
plt.xlim(0,40)
plt.axhline(0, linewidth=1, color ="k")
savefig("Earthquake_Record_vel_TZ_S2.png")
show()

plt.figure()
plt.title("Acceleration System 2 Tarzana")
plt.plot(tt,lista_u2p ,linewidth = 1.5,label="Newmark",alpha=0.7)
plt.plot(tt,u_2pi ,linewidth = 1.5,label="Duhamel's",alpha=0.7,color="g",linestyle="dashed")
plt.plot(tt,lista_u2p_cd,linestyle="dashed",linewidth = 1.5,label="Central difference",alpha=0.7)
plt.legend(loc="upper right",fontsize="larger")
plt.grid(b=True, alpha=0.7,linestyle="dashed", color="dimgray",linewidth=0.5)
plt.xlabel("Time [s]",size = 12)
plt.ylabel("Acceleration [m/s2]",size = 12)
plt.xlim(0,40)
plt.axhline(0, linewidth=1, color ="k")
savefig("Earthquake_Record_Acce_TZ_S2.png")
show()

###############################################################################

# /////////////////////////////////////////////////////////////////////////////


