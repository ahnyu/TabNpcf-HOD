import numpy as np
import math
from scipy import special
import glob
from numba import jit
from abacusnbody.data.compaso_halo_catalog import CompaSOHaloCatalog

###public HOD weight and Ngal
def weight(log_mhalo,hod):
    N_cent=np.zeros(log_mhalo.size)
    N_sate=np.zeros(log_mhalo.size)        
    log_mcut=hod[0]
    log_m1=hod[1]
    sigma=hod[2]
    alpha=hod[3]
    kappa=hod[4]
    m_halo=10.0**log_mhalo
    m_cut=10.0**log_mcut
    m_1=10.0**log_m1
    N_cent=np.array(1.0/2.0*special.erfc(np.log(m_cut/m_halo)/np.sqrt(2.)/sigma))
    mask=(m_halo>kappa*m_cut)
    N_sate[mask]=((m_halo[mask]-kappa*m_cut)/m_1)**alpha
    return N_cent,N_sate

def num_gal(hod,massbin):
    massbin_mid=massbin[2]
    num_massbin_halo=massbin[0]
    w=weight(massbin_mid,hod)
    return np.sum(w*num_massbin_halo)

def rsd_correction(cosmo):
    a=1/(1+cosmo['redshift'])
    om=cosmo['Om']
    ol=1.0-om
    return a,100.0*np.sqrt(om/a**3+ol)


def hal_cat(masscut,cosmo,pc_common,simulationIn,seed):
    allcat=sorted(glob.glob(simulationIn+'*.asdf'))
    ncat=len(allcat)
    hxall=[]
    hyall=[]
    hzall=[]
    pxall=[]
    pyall=[]
    pzall=[]
    nall=[]
    pindexall=[]
    pnumall=[]
    randhall=[]
    randpall=[]
    tmprsd=rsd_correction(cosmo)
    a=tmprsd[0]
    rsd=tmprsd[1]
    boxsize=pc_common['boxsize'][0]
    np.random.seed(seed)
    for i in range(ncat):
#        cat=CompaSOHaloCatalog(allcat[i],subsamples=dict(A=True,pos=True,vel=True))
        cat=CompaSOHaloCatalog(allcat[i],load_subsamples='A_halo_rv')

        massmask=(cat.halos['N']>=masscut)
        nhalos=len(cat.halos)
        nparts=len(cat.subsamples)
        rand_h=np.random.rand(nhalos)
        randhall.append(rand_h[massmask])
        rand_p=np.random.rand(nparts)
        randpall.append(rand_p)
        px=cat.subsamples['pos'][:,0]+1000.0
        py=cat.subsamples['pos'][:,1]+1000.0
        pz=cat.subsamples['pos'][:,2]+1000.0
        pvz=cat.subsamples['vel'][:,2]
        pznew=np.mod(pz+pvz/rsd/a,boxsize)
        pxall.append(px)
        pyall.append(py)
        pzall.append(pznew)
        hx=cat.halos['x_com'][:,0][massmask]+1000.0
        hy=cat.halos['x_com'][:,1][massmask]+1000.0
        hz=cat.halos['x_com'][:,2][massmask]+1000.0
        hvz=cat.halos['v_com'][:,2][massmask]
        hznew=np.mod(hz+hvz/rsd/a,boxsize)
        hxall.append(hx)
        hyall.append(hy)
        hzall.append(hznew)
        n=cat.halos['N'][massmask]
        nall.append(n)
        pindex=cat.halos['npstartA'][massmask]
        pindexall.append(pindex)
        pnum=cat.halos['npoutA'][massmask]
        pnumall.append(pnum)
    return hxall,hyall,hzall,pxall,pyall,pzall,nall,pindexall,pnumall,randhall,randpall


@jit(nopython=True, parallel=True)
def gal_cat_numba(hod,cosmo,halocat):
    hxall=halocat[0]
    hyall=halocat[1]
    hzall=halocat[2]
    pxall=halocat[3]
    pyall=halocat[4]
    pzall=halocat[5]
    nall=halocat[6]
    pindexall=halocat[7]
    pnumall=halocat[8]
    randhall=halocat[9]
    randpall=halocat[10]
    gx=[]
    gy=[]
    gz=[]
    log_mcut=hod[0]
    log_m1=hod[1]
    sigma=hod[2]
    alpha=hod[3]
    kappa=hod[4]
    mpart=cosmo['particleMass'][0]
    m_cut=10.0**log_mcut
    m_1=10.0**log_m1
    count_cent=0
    count_sate=0
    count_total=0
    for i in range(9):
        hmass=nall[i]*mpart
        nhalos_cat=len(nall[i])
        nparts_cat=len(pxall[i])
        for j in range(nhalos_cat):
            N_cent=1.0/2.0*math.erfc(np.log(m_cut/hmass[j])/sigma/np.sqrt(2.0))
            if(randhall[i][j] < N_cent):
                gx.append(hxall[i][j])
                gy.append(hyall[i][j])
                gz.append(hzall[i][j])
                count_cent+=1
                count_total+=1
            if(hmass[j]>kappa*m_cut):
                N_sate=((hmass[j]-kappa*m_cut)/m_1)**alpha
                nparts=pnumall[i][j]
                if(nparts>0):
                    sp=N_sate/nparts
                    for k in range(nparts):
                        p_ind=int(pindexall[i][j]+k)
                        if(randpall[i][p_ind] < sp):
                            count_sate+=1
                            count_total+=1
                            gx.append(pxall[i][p_ind])
                            gy.append(pyall[i][p_ind])
                            gz.append(pzall[i][p_ind])
    agx=np.array(gx)
    agy=np.array(gy)
    agz=np.array(gz)

    return agx,agy,agz


@jit(nopython=True)
def get_particle_catalog(halocat):
    pxall=halocat[3]
    pyall=halocat[4]
    pzall=halocat[5]
    nall=halocat[6]
    pindexall=halocat[7]
    pnumall=halocat[8]
    pxout=[]
    pyout=[]
    pzout=[]
    pnout=[]
    for ipart in range(9):
        px=pxall[ipart]
        py=pyall[ipart]
        pz=pzall[ipart]
        n=nall[ipart]
        pindex=pindexall[ipart]
        pnum=pnumall[ipart]
        nhalo=len(n)
        for ihalo in range(nhalo):
            if(pnum[ihalo]!=0):
                for i in range(pnum[ihalo]):
                    pxout.append(px[pindex[ihalo]+i])
                    pyout.append(py[pindex[ihalo]+i])
                    pzout.append(pz[pindex[ihalo]+i])
                    pnout.append(n[ihalo])
    pxout=np.array(pxout)
    pyout=np.array(pyout)
    pzout=np.array(pzout)
    pnout=np.array(pnout)
    return pxout,pyout,pzout,pnout


