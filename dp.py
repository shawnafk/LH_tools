import numpy as np
import matplotlib.pyplot as plt
import math
#version with out normalization
#fundamental plasma parameters
#constants
#Electron and proton charge: coluomb as unit_q
qe=1.6022e-19
qp=1.6022e-19
#electron and proton mass: kg as unit_m
me=9.1094e-31 
mp=1836*me
#permittivity in vaccum,F/m
epsi0=8.8543e-12
#permeability in vaccum H/m
mu0=4 * np.pi * 1e-7
c=3e8

#input parameters
#magnetic induction
b0 = 2
#number density
ne = 1.5e19 
ni = 1.5e19
#ion mass and charge ratio to a proton
Rm=2; Rq=1;
#qi and mi are:
mi=Rm*mp
qi=Rq*qp
#electron gyro and plasma frequency
wce = qe*b0/me
wpe = ((ne*qe**2)/(epsi0*me))**0.5
#ion gyro and plasma frequency
wci = qi*b0/mi
wpi = ((ni*qi**2)/(epsi0*mi))**0.5

#lower hybrid wave dispersion relation under omega_ci << w << omega_ce
#direct tensor epsi 
#coeffiency
#m-1
X=[]
Y=[]
#wlh=omega_pi/(1+omega_pe**2/omega_ce**2)**0.5
#k_para=2*np.pi/(1024*5e-4)
k_para=2*np.pi/(192*2e-4)
OMEGA=np.linspace(0.011,310.4,2000)*1e8
#solve ???
for w in OMEGA:
    R=1 - (wpe**2/w/(w-wce) + wpi**2/w/(w+wci))
    L=1 - (wpe**2/w/(w+wce) + wpi**2/w/(w-wci))
    P=1 - (wpe**2/w**2 + wpi**2/w**2)
    S = 0.5 * (R + L)
    D = 0.5 * (R - L)
    A=S
    B=R*L + P*S - P * (k_para*c/w)**2 - S * (k_para*c/w)**2
    C=P * (R*L - 2 * S * (k_para*c/w)**2 + (k_para*c/w)**4)
    #P6=-( 3 * (omega_pi*vti/omega)**2 + 3./4. * (omega_pe*omega*vte/omega_ce/omega_ce)**2)
    #print(P6)
    for _ in np.roots([A,-B,C]):
        if np.isreal(_):
        	X.append((np.real(_)))
	        Y.append(w)
       #	 X.append(w)
#print(X)
plt.plot(X,Y,'.',markersize=1.5)
#plt.plot(X,Y)
plt.show()
