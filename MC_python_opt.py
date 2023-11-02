 #####################################################################
#                                                                    #
#       MONTE CARLO NVT-CONSTANTES PARA ATOMOS LENNARD JONES         #
#                                                                    #
#                                                                    #
######################################################################
#  PROGRAMA DE SIMULACION DE MONTECARLO EN UN ENSAMBLE NVT-CONSTANTES
#
#  ESTE PROGRAMA TOMA UNA CONFIGURACON DE ATOMOS DE LENNARD JONES
#
#  Y REALIZA UNA SIMULACION CONVENCIONAL DE MC NVT. LA CAJA ES DE
#
#  LONGITUD UNITARIA, -0.5 A 0.5 Y NO HAY TABLAS DE CONSULTA
#
#
#  VARIABLES PRINCIPALES:
#
#
#  ENTERO   n                       NUMERO DE MOLECULAS
#
#  ENTERO   nstep                   MAXIMO NUMERO DE CICLOS
#
#  REAL     rx[n],ry[n],rz[n]       POSICIONES
#
#  REAL     dens                    DENSIDAD REDUCIDA
#
#  REAL     temp                    TEMPERATURA REDUCIDA
#
#  REAL     sigma                   DIAMETRO LJ REDUCIDO
#
#  REAL     rmin                    SEPARACION DE PARES MINIMA REDUCIDA
#
#  REAL     rcut                    RADIO DE CORTE REDUCIDO
#
#  REAL     drmax                   DESPLAZAMIENTO MAXIMO REDUCIDO
#
#  REAL     v                       ENERGIA POTENCIAL
#
#  REAL     w                       VIRIAL
#
#  REAL     press                   PRESION
#
#
# USO:
#
#
# EL PROGRAMA USA UNIDADES LENNARD JONES PARA LAS ENTRADAS DEL USARIO Y
#
# LAS SALIDAS PERO LLEVA A CABO LA SIMULACION EN UNA CAJA DE LONGITUD
#
# UNITARIA, POR EJEMPLO, PARA UNA LONGITUD DE CAJA L, Y PARAMETROS DE
#
# LENNARD-JONES EPSILUM Y SIGMA, LAS UNIDADES SON:
#
#  PROPIEDAD        UNIDADES LJ             UNIDADES DEL PROGRAMA
#
#  TEMPERATURA      EPSILUM/K               EPSILUM/K
#
#  PRESION          EPSILUM/SIGMA**3        EPSILUM
#
#  V                EPSILUM                 EPSILUM
#
#  DENS             1/SIGMA**3              1/L**3
#
#
# REFERNCIA DE LAS FUNCIONES:
#
#
# def sumup():
#
#       Calcula la energía potencial total para una configuración
#
# def energy():
#
#       Calcula la energía potencial del átomo i con todos los otros
#
#       átomos en el líquido
#
# def readcn():
#
#       Lee la configuración inicial
#
# def writcn():
#
#       Escribe la configuración de salida
#
#######################################################################

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

        rxij = rx[i] - rx[i+1:n+1]
        ryij = ry[i] - ry[i+1:n+1]
        rzij = rz[i] - rz[i+1:n+1]

#       mínima imagen de la separación de pares

        rxij = rxij -  np.round(rxij*boxinv)*box
        ryij = ryij -  np.round(ryij*boxinv)*box
        rzij = rzij -  np.round(rzij*boxinv)*box
        rijsq = rxij**2+ryij**2+rzij**2

        superp = len(rijsq[rijsq<rmin])

        if superp>0:

            ovrlap = True
            return(ovrlap)

        sr2 = sigsq/rijsq[rijsq<rcutsq]
        sr6 = sr2**3
        vij = sr6*(sr6-1)
        wij = sr6*(sr6-0.5)
        vf = vf+np.sum(vij)
        wf = wf+np.sum(wij)

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
    rxg = np.delete(rx,[0,ip])
    ryg = np.delete(ry,[0,ip])
    rzg = np.delete(rz,[0,ip])

#   bucle sobre todas las moléculas excepto i

    rxij = rxif-rxg
    ryij = ryif-ryg
    rzij = rzif-rzg

    rxij = rxij-np.round(rxij*boxinv)*box
    ryij = ryij-np.round(ryij*boxinv)*box
    rzij = rzij-np.round(rzij*boxinv)*box

    rijsq = rxij**2 + ryij**2 + rzij**2

    sr2 = sigsq/rijsq[rijsq<rcutsq]
    sr6 = sr2**3
    vij = sr6*(sr6-1)
    wij = sr6*(sr6-0.5)
    vf = np.sum(vij)
    wf = np.sum(wij)
    vf = 4*vf
    wf = 48*wf/3
    return(vf,wf)


#######################################################################
#       BUCLES SOBRE TODOS LOS CICLOS Y TODAS LAS MOLECULAS
#
#######################################################################

for istep in range(1,nstep+1):

    for i in range(1,n+1):

        rxiold = rx[i]
        ryiold = ry[i]
        rziold = rz[i]

        vold,wold=energy(rxiold,ryiold,rziold,i)

#       mover i y recoger la imagen central

        rxinew = rxiold + (2*np.random.rand()-1)*drmax
        ryinew = ryiold + (2*np.random.rand()-1)*drmax
        rzinew = rziold + (2*np.random.rand()-1)*drmax

        rxinew = rxinew - round(rxinew*boxinv)*box
        ryinew = ryinew - round(ryinew*boxinv)*box
        rzinew = rzinew - round(rzinew*boxinv)*box

#       calcular la energía de i en la nueva configuración

        vnew,wnew=energy(rxinew,ryinew,rzinew,i)

#       Comprobar la aceptación

        deltv = vnew - vold
        deltw = wnew - wold
        deltvb = beta * deltv

        if(deltv <= 0):

            v = v+deltv
            w = w+deltw
            rx[i] = rxinew
            ry[i] = ryinew
            rz[i] = rzinew
            acatma = acatma+1

        elif ( math.exp(-deltvb) > np.random.rand() ):

            v = v+deltv
            w = w+deltw
            rx[i] = rxinew
            ry[i] = ryinew
            rz[i] = rzinew
            acatma = acatma+1

        acm = acm+1

#       calcular los valores instantaneos

        vn = v/n
        press = dens*temp+w/vol
        press = press*sigma**3

#       promedios acomulados

        acv = acv+vn
        acp = acp+press
        acvsq = acvsq+vn**2
        acpsq = acpsq+press**2

###################################################################
#        FIN DEL CICLO SOBRE TODOS LOS ATOMOS
#
###################################################################
    if( istep%iratio == 0 ):

        ratio = acatma/(n*iratio)

        if(ratio >0.5):

            drmax = drmax*1.05

        else:

            drmax = drmax*0.95

        acatma=0

#Se calcula la N(r+Dr)
    if istep%ngr == 0:
        box2 = 0.5*box

        for i in range(1,n):

            dx = rx[i] - rx[i+1:n+1]
            dy = ry[i] - ry[i+1:n+1]
            dz = rz[i] - rz[i+1:n+1]

            dx = dx - np.round(dx*boxinv)*box
            dy = dy - np.round(dy*boxinv)*box
            dz = dz - np.round(dz*boxinv)*box

            rij = np.sqrt(dx**2+dy**2+dz**2)

            irij = np.round(rij[rij<=box2]/deltar)
            nga[irij.astype(int)]= nga[irij.astype(int)] + 2

# Escrbir la información de salida
    if (istep%iprint == 0):
    #if (istep%iratio == 0):
        print("           "+str(int(acm))+"     "+str(round(ratio,4))+"    "+str(round(vn,4))+"    "+str(round(press,4)))

####################################################################
#       FIN DEL BUCLE SOBRE TODOS LOS CICLOS
#
####################################################################

print("FIN DE LA CADENAD DE MARKOV ")

# Se calcula la rdf

nstepg = nstep/ngr
file2=open("rdf.dat","w")
if nstepg >0:

    nig = int(0.5*box/deltar)
    vol = box**3
    rhoj = n/vol
    denx = float(nstepg*n)
    fact = 4*np.pi*rhoj/3

    for jj in range(1,nig+1):
        r0[jj] = jj*deltar
        r = r0[jj] + deltar
        den = fact*(r**3-r0[jj]**3)
        gar = int(nga[jj])/den
        gr[jj] = gar/denx
        file2.write("{r0w:<25.20f} {grjj:<25.20f}\n".format(r0w=r0[jj],grjj=gr[jj]))
file2.close()

# graficar la rdf

xr=np.delete(r0,range(nig,10000))
yg=np.delete(gr,range(nig,10000))

plt.title("Función de distribución radial")
plt.xlabel("r")
plt.ylabel("g(r)")
plt.plot(xr,yg)
plt.savefig("rdf.png")
plt.show()

# CALCULAR Y ESCRIBIR LOS PROMEDIOS DE LA CORRIDA

avv = acv/acm
acvsq = (acvsq/acm)-avv**2
avp = acp/acm
acpsq = (acpsq/acm)-avp**2

print("Promedios")
print("<V/N>= ",round(avv,4))
print("<P>= ",round(avp,4))

# CALCULAR LAS FLUCTUACIONES

if acvsq>0:
    flv=acvsq**0.5
if acpsq>0:
    flp=acpsq**0.5

print("Fluctuaciones")
print("Fluctuación en <V/N> = ",round(flv,4))
print("Fluctuaciones en <P> = ",round(flp,4))
print("Fin de la simulación")
