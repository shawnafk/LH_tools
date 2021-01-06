import numpy as np
import matplotlib.pyplot as plt
from sympic_io.tools import unit as unit
class LH_Dp:
    def __init__(self,n,B,m,q,U):
        self.n=n
        self.B=B
        self.U=U
        self.m=m
        self.q=q
        self.w_ce =  unit.Oce(B,U)
        self.w_pe =  unit.Ope(n,U)
        self.w_pi  = unit.Opi(m,q,n,U)
        self.w_ci  = unit.Oci(m,q,B,U)
        self.w_gmg=np.sqrt(self.w_ci*self.w_ce)
        #ref this works under wci much less than wpi, i.e. high density low magnetic field
        self.wlh=self.w_pi/(1+self.w_pe**2/self.w_ce**2)**0.5
        #full  exp
        #wlh=(omega_pi**2/(1+omega_pe**2/omega_ce**2) + omega_ci**2)**0.5
    def omegakp(self,k_para,omega):
        w_ce=self.w_ce
        w_ci=self.w_ci
        w_pe=self.w_pe
        w_pi=self.w_pi
        kv=[]
        c=1
        n_para=k_para*c/omega
        epsi_perp=1+(w_pe/w_ce)**2 - (w_pi/omega)**2
        epsi_para=1-(w_pe/omega)**2 - (w_pi/omega)**2
        epsi_xy=(w_pe)**2 / (w_ce*omega)
        P0=epsi_para*( (n_para**2-epsi_perp)**2 - epsi_xy**2)
        P2=(epsi_perp + epsi_para)*(n_para**2 - epsi_perp) + epsi_xy**2
        P4=epsi_perp
        for _ in np.roots([P4,P2,P0]):
            if(np.isreal(_) and _ > 0):
                kv.append(np.sqrt(_)*omega/c)
        #return n*w/c = k
        return kv

    #lower hybrid wave dispersion relation under omega_ci << w << omega_ce
    #direct tensor epsi 
    #coeffiency
    #m-1
    def omegakp_hot(self,k_para,omega,vte):
        c=1
        w_ce=self.w_ce
        w_ci=self.w_ci
        w_pe=self.w_pe
        w_pi=self.w_pi
        kv=[]
        n_para=k_para*c/omega
        epsi_perp=1+(w_pe/w_ce)**2 - (w_pi/omega)**2
        epsi_para=1-(w_pe/omega)**2 - (w_pi/omega)**2
        epsi_xy=(w_pe)**2 / (w_ce*omega)
        P0=epsi_para*( (n_para**2-epsi_perp)**2 - epsi_xy**2)
        P2=(epsi_perp + epsi_para)*(n_para**2 - epsi_perp) + epsi_xy**2
        P4=epsi_perp
        P6=0.75*(w_pe*omega/w_ce**2)**2*vte**2
        for _ in np.roots([P6,P4,P2,P0]):
            if(np.isreal(_) and _ > 0):
                kv.append(np.sqrt(_)*omega/c)
        return kv
    #static k vert >> k para
    #warm simple
    def omegakp_esh(self,k_para,omega,vte):
        kv=[]
        P0=self.m*k_para**2
        P2=1-omega**2/wlh**2
        P4=0.75*vthe**2/omega_ce**2
        for _ in np.roots([P4,P2,P0]):
            if(np.isreal(_) and _ > 0):
                kv.append(np.sqrt(_-k_para**2))
        return kv
    def omegakp_es(self,k_para,omega,vte):
        kv=[]
        kperp_pre = kL(k_para,omega,omega_pe,omega_ce,omega_pi)
        if(np.isreal(kperp_pre) and kperp_pre > 0):
            kv.append(np.sqrt(kperp_pre))
        return kv

    #from plot.py, general functions
    #kvert >> kpara 
    def kvert_es(self,kpp,omega):
        w_ce=self.w_ce
        w_ci=self.w_ci
        w_pe=self.w_pe
        w_pi=self.w_pi
        kvt=1-(w_pi/omega)**2+(w_pe/w_ce)**2
        k11=1-(w_pe/omega)**2
        return np.sqrt(-(kpp**2*k11/kvt))

    def gamma(self,kpara,kvert,omega,vt):
        w_ce=self.w_ce
        w_ci=self.w_ci
        w_pe=self.w_pe
        w_pi=self.w_pi
        K=(kpara**2 + kvert**2)**0.5
        ld=vt/wpe
        gamma=sqrt(pi)/K**2/ld**2*omega/2**0.5/kpara/vt * exp(-(omega/kpara/vt)**2/2)/(kvert**2/K**2 * w_pi**2/omega**3 + kpara**2/K**2*w_pe**2/omega**3)/2
        return gamma
    def wz(self,kp,omega):
        w_ce=self.w_ce
        w_ci=self.w_ci
        w_pe=self.w_pe
        w_pi=self.w_pi
        ne=self.n
        ni=ne/self.q
        kxx=1-(wpi/w)**2 + (wpe/wce)**2
        kzz=1-(wpe/w)**2;
        return 2*np.pi/(kp*np.sqrt(abs(kzz/kxx)))
    #opt density
    def opt_n(self,kp):
        U=self.U
        q=1.6e-19/U[0]['CHARGE']
        return kp**2/q**2 / U[0]['LENGTH']**3
    #critical refaction index
    def crt_n(self):
        w_ce=self.w_ce
        w_ci=self.w_ci
        w_pe=self.w_pe
        w_pi=self.w_pi
        w_gmg=self.w_gmg
        nc=sqrt(1+(w_pi/omega)**2*((omega/w_gmg)**2-(omega**2/(omega**2-w_ci**2)))) + w_pi/w_gmg
        return nc
    def vg(self,kp,w):
        wce=self.w_ce
        wci=self.w_ci
        wpe=self.w_pe
        wpi=self.w_pi
        return (-1.*(-((kp**2*(1 - wpe**2/w**2))/(1 + wpe**2/wce**2 - wpi**2/w**2)))**0.5* (1.*w**2*wce**2 + 1.*w**2*wpe**2 - 1.*wce**2*wpi**2)**2)/ (kp**2*w*wce**2*(1.*wce**2*wpe**       2 + 1.*wpe**4 - 1.*wce**2*wpi**2))

'''
    def gamma_col(self,11,kL,wr,nue,nui,wpe,wpi,wce):
    K=(k11**2 + kL**2)**0.5
        gamma=( ( (wpe/wr)**2 * (k11/K)**2 + (wpe/wce)**2 * (kL/K)**2 ) *nue + (wpi/wr)**2*nui)/((wpe/wr)**2 * (k11/K)**2 + (wpi/wr)**2)
        return gamma
    #return vertical wave length with given karallel
'''

