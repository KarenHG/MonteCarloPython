import numpy as np
import sys
import math
import matplotlib.pyplot as plt

# Leer los datos de entrada

print("PROGRAMA MCLJ\n")
print("PROGRAMA MONTECARLO NVT CONSTANTES\n")
print("PARA ATOMOS DE TIPO LENNARD JONES\n")
title=input("Ingrese el título de la corrida\n")
nstep=int(input("Ingrese el número de ciclos\n"))
iprint=int(input("Ingrese el número de pasos entre las lineas de salida\n"))
ngr=int(input("Ingrese la frecuencia rdf\n"))
iratio=int(input("Ingrese el intervalo para actualizar el desplazamiento máximo\n"))
cnfile=input("Ingrese el nombre del documento de la configuración inicial\n")
print("Ingrese lo siguiente en unidades Lennard-Jones\n")
dens=float(input("Ingrese la densidad\n"))
temp=float(input("Ingrese la temperatura\n"))
rcut=float(input("Ingrese el rádio de corte\n"))

# Leemos la configuración inicial

file = open(cnfile,"r")
a=file.readline().split()
n=int(a[0])
box=float(a[1])

rx=np.zeros(n+1)
ry=np.zeros(n+1)
rz=np.zeros(n+1)

for i in range(1,n+1):
    partes = file.readline().split()
    rx[i]=float(partes[0])
    ry[i]=float(partes[1])
    rz[i]=float(partes[2])
file.close()

boxinv=1/box

#escribir los datos de entrada

print("\n ",title)
print(" Número de átomos                                 " ,n)
print(" Número de ciclos                                 " ,nstep)
print(" Frecuecia de salida                              " ,iprint)
print(" Frecuencia de actualizacion de radio             ",iratio)
print(" Nombre del documento de la configuración inicial ",cnfile)
print(" Temperatura                                      ",temp)
print(" Densidad                                         ",dens)
print(" Potencial de corte                               ",rcut)

vol=box**3
rho=n/vol

print(" Volumen                                          ",vol)
print(" Rho                                              ",rho)
print(" boxx,boxy,boxz                                   ",round(box,3),round(box,3),round(box,3))

#convertir los datos de entrda a unidades del programa

beta=1/temp
sigma=1
rmin=0.70*sigma
rcut=rcut*sigma
drmax=0.15*sigma
denslj=dens
dens=dens/(sigma**3)

if rcut > 0.5*box:
    sys.exit("\n RADIO DE CORTE DEMASIADO LARGO")

# Acomuladores en cero

acv=0
acvsq=0
acp=0
acpsq=0
flv=0
flp=0
acm=0
acatma=0
deltar=0.05

#revisar valor de nmax
nga=np.zeros(10000)
gr=np.zeros(10000)
r0=np.zeros(10000)

# escribir información util

print(" sigma/box                                        ",sigma)
print(" rmin/box                                         ",rmin)
print(" rcut/box                                         ",rcut)

def sumup(rcutf,rminf,sigmaf):

##################################################################
# CALCULA LA ENERGIA POTENCIAL DE LA CONFIGURACION
#
# VARIABLES PRINCIPALES:
#
# ENTERO        n                   EL NUMERO DE ATOMOS
#
# REAL          rx[n],ry[n],rz[n]   LA POSICION DE LOS ATOMOS
#
# REAL          vf                  LA ENERGIA POTENCIAL
#
# REAL          wf                  EL VIRIAL
#
# BOLEANO       ovrlap              VARIABLE DE PRUEBA PARA
#                                   SUPERPOSICION DE LOS ATOMOS
#
# USO:
#
# LA FUNCION RETORNA LA ENERGIA POTENCIAL TOTAL AL INICIO Y AL
#
# FINAL DE LA CORRIDA
#
###################################################################

    ovrlapf = False
    rcutsq = rcutf**2
    rminsq = rminf**2
    sigsq = sigmaf**2

    vf=0
    wf=0

#   bucle sobre todos los pares in líquido

    for i in range(1,n):

        rxi = rx[i]
        ryi = ry[i]
        rzi = rz[i]

        for j in range(i+1,n+1):

            rxij = rxi - rx[j]
            ryij = ryi - ry[j]
            rzij = rzi - rz[j]

#           mínima imagen de la separación de pares

            rxij = rxij -  round(rxij*boxinv)*box
            ryij = ryij -  round(ryij*boxinv)*box
            rzij = rzij -  round(rzij*boxinv)*box
            rijsq = rxij**2+ryij**2+rzij**2

            if rijsq < rminsq:

                ovrlap = True
                return(ovrlap)

            elif rijsq < rcutsq:

                sr2 = sigsq/rijsq
                sr6 = sr2**3
                vij = sr6*(sr6-1)
                wij = sr6*(sr6-0.5)
                vf = vf+vij
                wf = wf+wij

    vf = 4*vf
    wf = 48*wf/3
    return(vf,wf)

#calcular la energía inicial y comprobar si hay superposicion

if sumup(rcut,rmin,sigma)==True:
    sys.exit("\n SOBREPOSICION EN LA CONFIGURACION INICIAL")

v,w=sumup(rcut,rmin,sigma)
vvv = v
vs = v/n
ws = w/vol
ps = dens*temp+w/vol
ps = ps*sigma**3

print(" V inicial                                        ",round(vs,4))
print(" W inicial                                        ",round(ws,4))
print(" P inicial                                        ",round(ps,4))
print("             INICIO DE LA CADENA DE MARKOV        ")
print("            NMOV     RATIO       V/N     P        ")

def energy(rxif,ryif,rzif,ip):
    rcutsq = rcut**2
    sigsq = sigma**2
    vf = 0
    wf = 0

#   bucle sobre todas las moléculas excepto i

    for j in range(1,n+1):

        if ip!=j:

            rxij = rxif-rx[j]
            ryij = ryif-ry[j]
            rzij = rzif-rz[j]

            rxij = rxij-round(rxij*boxinv)*box
            ryij = ryij-round(ryij*boxinv)*box
            rzij = rzij-round(rzij*boxinv)*box

            rijsq = rxij**2 + ryij**2 + rzij**2

            if rijsq < rcutsq:

                sr2 = sigsq/rijsq
                sr6 = sr2**3
                vij = sr6*(sr6-1)
                wij = sr6*(sr6-0.5)
                vf = vf+vij
                wf = wf+wij
    vf = 4*vf
    wf = 48*wf/3
    return(vf,wf)


#######################################################################
#       BUCLES SOBRE TODOS LOS CICLOS Y TODAS LAS MOLECULAS
#
#######################################################################

#for istep in range(1,nstep+1):
for istep in range(1,3):

    #for i in range(1,n+1):
    for i in range(1,3):

        rxiold = rx[i]
        ryiold = ry[i]
        rziold = rz[i]

        vold,wold=energy(rxiold,ryiold,rziold,i)
        print(vold,wold)
