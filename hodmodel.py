import numpy as np
import math
from scipy import special
import glob
from numba import jit
from abacusnbody.data.compaso_halo_catalog import CompaSOHaloCatalog

#weight
def weight_lrg(log_mhalo,params_lrg,model_lrg):      
    
    
    m_halo=10.0**log_mhalo
    if(model_lrg['cent']==1):
        log_mcut=params_lrg['logMcut']
        sigma=params_lrg['sigma']
        pmax=params_lrg['pmax']
        m_cut=10.0**log_mcut
        N_cent=Ncent(log_mhalo,log_mcut,sigma,pmax) 
    if(model_lrg['sate']==1):
        log_m1=params_lrg['logM1']
        alpha=params_lrg['alpha']
        kappa=params_lrg['kappa']
        m_1=10.0**log_m1
        N_sate=Nsate(m_halo,m_cut,m_1,alpha,kappa)
    return N_cent,N_sate

def weight_elg(log_mhalo,params_elg,model_elg):
    
    
    m_halo=10.0**log_mhalo
    if(model_elg['cent']==1):
        log_mcut=params_elg['logMcut']
        sigma=params_elg['sigma']
        pmax=params_elg['pmax']
        m_cut=10.0**log_mcut
        N_cent=Ncent(log_mhalo,log_mcut,sigma,pmax)
    elif(model_elg['cent']==2):
        log_mcut=params_elg['logMcut']
        sigma=params_elg['sigma']
        pmax=params_elg['pmax']
        gamma=params_elg['gamma']
        Q=params_elg['Q']
        maxpdf=params_elg['maxpdf']
        m_cut=10.0**log_mcut
        N_cent=Ncent_elg_hmq(log_mhalo,log_mcut,sigma,pmax,gamma,Q,maxpdf)

    if(model_elg['sate']==1):
        log_m1=params_elg['logM1']
        alpha=params_elg['alpha']
        kappa=params_elg['kappa']
        m_1=10.0**log_m1
        N_sate=Nsate(m_halo,m_cut,m_1,alpha,kappa)
    return N_cent,N_sate

def weight_qso(log_mhalo,params_qso,model_qso):
    

    m_halo=10.0**log_mhalo
    if(model_qso['cent']==1):
        log_mcut=params_qso['logMcut']
        sigma=params_qso['sigma']
        pmax=params_qso['pmax']
        m_cut=10.0**log_mcut
        N_cent=Ncent(log_mhalo,log_mcut,sigma,pmax)

    if(model_qso['sate']==1):
        log_m1=params_qso['logM1']
        alpha=params_qso['alpha']
        kappa=params_qso['kappa']
        m_1=10.0**log_m1
        N_sate=Nsate(m_halo,m_cut,m_1,alpha,kappa)
    elif(model_qso['sate']==2):
        log_m1=params_qso['logM1']
        alpha=params_qso['alpha']
        log_mmin=params_qso['logMmin']
        m_1=10.0**log_m1
        m_min=10.0**log_mmin
        N_sate=Nsate_qso(m_halo,m_cut,m_1,alpha,m_min)
    return N_cent,N_sate

#Ncent Nsate

def Ncent(log_mh,log_mc,sigma,pmax):   
    return pmax*1.0/2.0*special.erfc((log_mc-log_mh)/np.sqrt(2)/np.log10(np.e)/sigma)

def Ncent_elg_hmq(log_mh,log_mc,sigma,pmax,gamma,Q,maxpdf):
    A=(pmax-1./Q)/(maxpdf)
    stepfunc=1./2./Q*(1+special.erf((log_mh-log_mc)/0.01))
    return A*pdf_normal(log_mh,log_mc,sigma,gamma)+stepfunc

def pdf_normal(x,mu,sigma,alpha):
    normal=np.exp(-1./2.*np.power((x-mu)/sigma,2))*(1/sigma/np.sqrt(2.*np.pi))
    cdf=1./2.*(1+special.erf((alpha*(x-mu))/np.sqrt(2.)/sigma))
    pdf=2*normal*cdf
    return pdf

def Nsate(Mh,Mc,M1,alpha,kappa):
    N_sate=np.zeros(len(Mh))
    mask=(Mh>kappa*Mc)
    N_sate[mask]=((Mh[mask]-kappa*Mc)/M1)**alpha
    return N_sate

def Nsate_qso(Mh,Mc,M1,alpha,Mmin):
    return (Mh/M1)**alpha*np.exp(-Mmin/Mh)

def findmaxpdf(log_mc,sigma,gamma):
    tmp_logmh=np.linspace(10.,15.5,1000)
    tmp_pdf=[]
    for i in range(1000):
        tmp_pdf.append(pdf_normal(tmp_logmh[i],log_mc,sigma,gamma))
    return np.max(np.array(tmp_pdf))


def num_elg(hod,bin_mid,num_halo,model_elg):
    w=weight_elg(bin_mid,hod,model_elg)
    return np.sum(w*num_halo)
def num_qso(hod,bin_mid,num_halo,model_qso):
    w=weight_qso(bin_mid,hod,model_qso)
    return np.sum(w*num_halo)
def num_lrg(hod,bin_mid,num_halo,model_lrg):
    w=weight_lrg(bin_mid,hod,model_lrg)
    return np.sum((w[0]+w[0]*w[1])*num_halo)
