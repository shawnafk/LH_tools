#!/usr/bin/env python3
# coding=utf-8
import numpy as np
ne=np.linspace(0.1,1,100)*1e14
ni=np.linspace(0.1,1,100)*1e14
mu=3672
B=20000
z=1
wpe=5.64e4*ne**0.5
wpi=1.32e3*z*mu**-0.5*ni**0.5
wce=1.76e7*B
w=4*np.pi*2.45e9
k_para=2*np.pi/(192*2e-2)
k_perp=(-((1-(wpe/w)**2)/(1-(wpi/w)**2 + (wpe/wce)**2) * k_para**2))**0.5

k_vert=np.sqrt(((omega*omega_pe**2/omega_ce/c**2)**2-k_para**4)/k_para**2)
