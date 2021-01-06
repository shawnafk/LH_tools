#!/usr/bin/env python3
# coding=utf-8
#Di=
#De=
#b90=
#Coulomb_log=
import numpy as np
import matplotlib.pyplot as plt
CL=15
mu=0.5
n=1e13
Te=np.linspace(10,15,1000)*1e3
Ti=np.linspace(1,10,1000)*1e3
nu_ee=2.9e-6*n*CL*Te**-1.5
nu_ii=4.8e-8*n*CL*Ti**-1.5*mu**-0.5

#vf=np.linspace(0.2,0.35,100)
vf=0.513
e_f=0.5*511e3*vf**2
#vs=np.linspace(0.001,0.002,100)
vs=0.8/1024/2
e_s=0.5*511e3*vs**2
nu_ei_f=2.1e-9*mu**-1*Ti*e_f**-2.5 * n *CL
nu_ee_f=3.9e-6*Ti*e_f**-2.5 * n *CL

nu_ei_s=2.5e-4*mu**0.5*Ti**-0.5*e_s**-1 * n * CL
nu_ee_s=2.9e-6*Ti**-0.5*e_s**-1 * n * CL


plt.semilogy(Te,nu_ee_f)
plt.semilogy(Ti,nu_ee/2)
plt.legend(['nu_ee_fast','nu_ei'])
#plt.semilogy(Te,nu_ee_s)
#plt.semilogy(Ti,nu_ei_f)
#plt.semilogy(Ti,nu_ei_s,label='nu_ei_slow')
