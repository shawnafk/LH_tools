import numpy as np
import matplotlib.pyplot as plt
from lib_zjs import unit as unit
c=1
L=2e-4
n=1e18
B=0.5
U=unit.gen_unit(L)
omega_ce =  unit.Oce(B,U)
omega_pe =  unit.Ope(n,U)
omega_pi  = unit.Opi(1,1,n,U)
omega_ci  = unit.Oci(1,1,B,U)
#ref this works under wci much less than wpi, i.e. high density low magnetic field
wlh=omega_pi/(1+omega_pe**2/omega_ce**2)**0.5
#full  exp
wlh=(omega_pi**2/(1+omega_pe**2/omega_ce**2) + omega_ci**2)**0.5
#lower hybrid wave dispersion relation under omega_ci << w << omega_ce
#direct tensor epsi 
#coeffiency
#m-1
w=[]
kv=[]
#k_para=2*np.pi/(1024*5e-4)
#wid=9.5e-3
#n=wid*8/L
k_para=2*np.pi/(720)
OMEGA=np.linspace(0.000001,0.010,20000)
realomega=OMEGA*3e8/L
N_PARA=k_para*c/OMEGA
#solve ???
for n_para,omega in zip(N_PARA,OMEGA):
#[n_para,omega] = [k_para*c/OMEGA,OMEGA]
	epsi_perp=1+(omega_pe/omega_ce)**2 - (omega_pi/omega)**2
	epsi_para=1-(omega_pe/omega)**2 - (omega_pi/omega)**2
	epsi_xy=(omega_pe)**2 / (omega_ce*omega)
	P0=epsi_para*( (n_para**2-epsi_perp)**2 - epsi_xy**2)
	P2=(epsi_perp + epsi_para)*(n_para**2 - epsi_perp) + epsi_xy**2
	P4=epsi_perp
	for _ in np.roots([P4,P2,P0]):
	    if(np.isreal(_) and _ > 0):
	        w.append(omega)
	        kv.append(np.sqrt(_)*omega/c)
	 #       print(kv)
DP=zip(kv,w)
res=sorted(DP,key= lambda v:v[0])
kv,w=zip(*res)
plt.semilogx(kv,w)
'''
w=[]
kv=[]
vte=0.00/3
for n_para,omega in zip(N_PARA,OMEGA):
#[n_para,omega] = [k_para*c/OMEGA,OMEGA]
	epsi_perp=1+(omega_pe/omega_ce)**2 - (omega_pi/omega)**2
	epsi_para=1-(omega_pe/omega)**2 - (omega_pi/omega)**2
	epsi_xy=(omega_pe)**2 / (omega_ce*omega)
	P0=epsi_para*( (n_para**2-epsi_perp)**2 - epsi_xy**2)
	P2=(epsi_perp + epsi_para)*(n_para**2 - epsi_perp) + epsi_xy**2
	P4=epsi_perp
	P6=0.75*(omega_pe*omega/omega_ce**2)**2*vte**2
	for _ in np.roots([P6,P4,P2,P0]):
	    if(np.isreal(_) and _ > 0):
	        w.append(omega)
	        kv.append(np.sqrt(_)*omega/c)
	 #       print(kv)
DP=zip(kv,w)
res=sorted(DP,key= lambda v:v[0])
kv,w=zip(*res)
plt.semilogx(kv,w)
plt.plot(kv,w)
plt.show()

#static k vert >> k para
#warm simple
w=[]
kv=[]
vthe=0.07/3
for omega in OMEGA:
#[n_para,omega] = [k_para*c/OMEGA,OMEGA]
	mi=1826;me=1
	P0=mi/me*k_para**2
	P2=1-omega**2/wlh**2
	P4=0.75*vthe**2/omega_ce**2
	for _ in np.roots([P4,P2,P0]):
	    if(np.isreal(_) and _ > 0):
	        w.append(omega)
	        kv.append(np.sqrt(_-k_para**2))
DP=zip(kv,w)
res=sorted(DP,key= lambda v:v[0])
kv,w=zip(*res)
plt.semilogx(kv,w)


w=[]
kv=[]
def kL(kpp,w,wpe,wce,wpi):
	kvt=1-(wpi/w)**2+(wpe/wce)**2
	k11=1-(wpe/w)**2
	return -(kpp**2*k11/kvt)
for omega in OMEGA:
	kperp_pre = kL(k_para,omega,omega_pe,omega_ce,omega_pi)
	if(np.isreal(kperp_pre) and kperp_pre > 0):
		w.append(omega)
		kv.append(np.sqrt(kperp_pre))
DP=zip(kv,w)
res=sorted(DP,key= lambda v:v[0])
kv,w=zip(*res)
plt.semilogx(kv,w)
plt.show()


kw=np.load('kv_w_pulse.npy')
x,y=kw.shape
x=int(x/2)
y=int(y/2)
X=np.linspace(0,np.pi,x)
Y=np.linspace(0,np.pi/(0.5*200),y)
K,W=np.meshgrid(X,Y)
plt.contour(K,W,abs(np.fft.fftshift(kw)[x:2*x,y:2*y]).transpose(),levels=3)
'''

#from plot.py, general functions
def gamma(k11,kL,wr,vt,wpe,wpi):
	K=(k11**2 + kL**2)**0.5
	ld=vt/wpe
	gamma=sqrt(pi)/K**2/ld**2*wr/2**0.5/k11/vt * exp(-(wr/k11/vt)**2/2)/(kL**2/K**2 * wpi**2/wr**3 + k11**2/K**2*wpe**2/wr**3)/2
	return gamma
def gamma_col(k11,kL,wr,nue,nui,wpe,wpi,wce):
	K=(k11**2 + kL**2)**0.5
	gamma=( ( (wpe/wr)**2 * (k11/K)**2 + (wpe/wce)**2 * (kL/K)**2 ) *nue + (wpi/wr)**2*nui)/((wpe/wr)**2 * (k11/K)**2 + (wpi/wr)**2)
	return gamma
def kL(kpp,w,wpe,wce,wci):
	kvt=1-(wpi/w)**2+(wpe/wce)**2
	k11=1-(wpe/w)**2
	return sqrt(-(kpp**2*k11/kvt))
def vg(kp,w,wpe,wce,wpi):
	return (-1.*(-((kp**2*(1 - wpe**2/w**2))/(1 + wpe**2/wce**2 - wpi**2/w**2)))**0.5* (1.*w**2*wce**2 + 1.*w**2*wpe**2 - 1.*wce**2*wpi**2)**2)/ (kp**2*w*wce**2*(1.*wce**2*wpe**       2 + 1.*wpe**4 - 1.*wce**2*wpi**2))




Ogmg=sqrt(Oci(m,B,U)*Oce(B,U))
nc=sqrt(1+(Opi(m,n,U)/omega)**2*((omega/Ogmg)**2-(omega**2/(omega**2-Oci(m,B,   U)**2)))) + Opi(m,n,U)/Ogmg


import matplotlib as mpl
mpl.rcParams['legend.edgecolor']= '1'
mpl.rcParams['legend.framealpha']= '0'
mpl.rcParams['xtick.direction']= 'in'
mpl.rcParams['ytick.direction']= 'in'
mpl.rcParams['xtick.labelsize']=15
mpl.rcParams['ytick.labelsize']=15
#mpl.rc('text', usetex=True)

def wz(w,kp,ne,B,mu,Z,U):
	ni=ne/Z;
	wci=Oci(mu,Z,B,U)
	wce=Oce(B,U)
	wpi=Opi(mu,Z,ni,U)
	wpe=Ope(ne,U)
	kxx=1-(wpi/w)**2 + (wpe/wce)**2
	kzz=1-(wpe/w)**2;
	return 2*np.pi/(kp*np.sqrt(abs(kzz/kxx)))
w=0.0043
dx=1e-3;U=gen_unit(dx)
mu=1;Z=1
#ne_pork=1e18*np.linspace(1e-2,1,1000)
B0=0.23;
B=np.linspace(B0,2.5,100)
kp=2*pi/(384)
wl=wz(w,kp,1e18,B,mu,Z,U)
#wl2=wz(w,kp,2*ne_pork,B,mu,Z,U)

#opt density
def opt_n(kp,U):
	q=1.6e-19/U[0]['CHARGE']
	return kp**2/q**2 / U[0]['LENGTH']**3

#plot(ne_pork,wl)
plot(B,wl)
show()
