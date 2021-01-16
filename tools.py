import numpy as np
import matplotlib.pyplot as plt
from sympic_io.tools import unit as unit
class LH_Dp:
    def __init__(self,m,q,U):
        self.U=U
        self.m=m
        self.q=q
        #ref this works under wci much less than wpi, i.e. high density low magnetic field
        #full  exp
        #wlh=(omega_pi**2/(1+omega_pe**2/omega_ce**2) + omega_ci**2)**0.5
    def pomega(self,n,B):
        w_ce =  unit.Oce(B,self.U)
        w_pe =  unit.Ope(n,self.U)
        w_pi  = unit.Opi(self.m,self.q,n,self.U)
        w_ci  = unit.Oci(self.m,self.q,B,self.U)
        w_gmg=np.sqrt(w_ci*w_ce)
        w_lh=w_pi/(1+w_pe**2/w_ce**2)**0.5
        print('w_ce={:f}'.format(w_ce))
        print('w_pe={:f}'.format(w_pe))
        print('w_ci={:f}'.format(w_ci))
        print('w_pi={:f}'.format(w_pi))
        print('w_gmg={:f}'.format(w_gmg))
        print('w_lh={:f}'.format(w_lh))
    
    def omegakv(self,n,B,k_para,omega):
        w_ce =  unit.Oce(B,self.U)
        w_pe =  unit.Ope(n,self.U)
        w_pi  = unit.Opi(self.m,self.q,n,self.U)
        w_ci  = unit.Oci(self.m,self.q,B,self.U)
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
    def omegakv_hot(self,n,B,k_para,omega,vte):
        c=1
        w_ce =  unit.Oce(B,self.U)
        w_pe =  unit.Ope(n,self.U)
        w_pi  = unit.Opi(self.m,self.q,n,self.U)
        w_ci  = unit.Oci(self.m,self.q,B,self.U)
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
    def omegakv_esh(self,n,B,k_para,omega,vte):
        w_ce =  unit.Oce(B,self.U)
        w_pe =  unit.Ope(n,self.U)
        w_lh=w_pi/(1+w_pe**2/w_ce**2)**0.5

        kv=[]
        P0=self.m*k_para**2
        P2=1-omega**2/w_lh**2
        P4=0.75*vte**2/w_ce**2
        for _ in np.roots([P4,P2,P0]):
            if(np.isreal(_) and _ > 0):
                kv.append(np.sqrt(_-k_para**2))
        return kv
    def omegakv_es(self,n,B,k_para,omega,vte):
        w_ce =  unit.Oce(V,self.U)
        w_pe =  unit.Ope(n,self.U)
        w_pi  = unit.Opi(self.m,self.q,n,self.U)
        kv=[]
        kperp_pre = kL(k_para,omega,w_pe,w_ce,w_pi)
        if(np.isreal(kperp_pre) and kperp_pre > 0):
            kv.append(np.sqrt(kperp_pre))
        return kv

    #from plot.py, general functions
    #kvert >> kpara 
    def kvert_es(self,n,B,kpp,omega):
        w_ce =  unit.Oce(B,self.U)
        w_pe =  unit.Ope(n,self.U)
        w_pi  = unit.Opi(self.m,self.q,n,self.U)
        kvt=1-(w_pi/omega)**2+(w_pe/w_ce)**2
        k11=1-(w_pe/omega)**2
        return np.sqrt(-(kpp**2*k11/kvt))

    def gamma(self,n,B,kpara,omega,vt):
        w_pe =  unit.Ope(n,self.U)
        w_pi  = unit.Opi(self.m,self.q,n,self.U)
        kvert= kvert_es(self,n,B,kpara,omega)
        K=(kpara**2 + kvert**2)**0.5
        ld=vt/wpe
        gamma=np.sqrt(np.pi)/K**2/ld**2*omega/2**0.5/kpara/vt * exp(-(omega/kpara/vt)**2/2)/(kvert**2/K**2 * w_pi**2/omega**3 + kpara**2/K**2*w_pe**2/omega**3)/2
        return gamma
#opt density
    def opt_n(self,kp):
        U=self.U
        q=1.6e-19/U[0]['CHARGE']
        return kp**2/q**2 / U[0]['LENGTH']**3
    #critical refaction index
#omega > w_gmg
    def crt_n(self,n,B,omega):
        w_ce =  unit.Oce(B,self.U)
        w_pi  = unit.Opi(self.m,self.q,n,self.U)
        w_ci  = unit.Oci(self.m,self.q,B,self.U)
        w_gmg=np.sqrt(w_ci*w_ce)
        nc=np.sqrt(1+(w_pi/omega)**2*((omega/w_gmg)**2-(omega**2/(omega**2-w_ci**2)))) + w_pi/w_gmg
        return nc
    def vg(self,n,B,kp,w):
        kvert= kvert_es(self,n,B,kp,w)
        w_ce =  unit.Oce(B,self.U)
        w_pe =  unit.Ope(n,self.U)
        w_pi  = unit.Opi(self.m,self.q,n,self.U)
        w_ci  = unit.Oci(self.m,self.q,B,self.U)
        vgv= (-1.*(-((kp**2*(1 - wpe**2/w**2))/(1 + wpe**2/wce**2 - wpi**2/w**2)))**0.5* (1.*w**2*wce**2 + 1.*w**2*wpe**2 - 1.*wce**2*wpi**2)**2)/ (kp**2*w*wce**2*(1.*wce**2*wpe**       2 + 1.*wpe**4 - 1.*wce**2*wpi**2))
        vgp= -kvert /kp * vgv
        return vgp,vgv

'''
    def gamma_col(self,11,kL,wr,nue,nui,wpe,wpi,wce):
    K=(k11**2 + kL**2)**0.5
        gamma=( ( (wpe/wr)**2 * (k11/K)**2 + (wpe/wce)**2 * (kL/K)**2 ) *nue + (wpi/wr)**2*nui)/((wpe/wr)**2 * (k11/K)**2 + (wpi/wr)**2)
        return gamma
    #return vertical wave length with given karallel
'''
'''
def wz(self,n,B,kp,omega):
        w_ce=self.w_ce
        w_ci=self.w_ci
        w_pe=self.w_pe
        w_pi=self.w_pi
        ne=self.n
        ni=ne/self.q
        kxx=1-(wpi/w)**2 + (wpe/wce)**2
        kzz=1-(wpe/w)**2;
        return 2*np.pi/(kp*np.sqrt(abs(kzz/kxx)))
'''

