import numpy as np
import ioconfig
#from ioconfig import dofitting,HODcatalog,params,wp2,xi3,fit_common,pc_common,HODcat,cosmo,paircountsIn,paircountsOut,triangleIn,simulationIn,chainsOut,covIn,data2In,data3In
from ioconfig import *
import hodmodel
from multiprocessing import Pool
import emcee
import paircounts as pc
import counts2clustering as c2c
import analyticalRandom as ar


###MCMC
def log_probability(theta):
    log_mcut,log_m1,sigma,alpha,kappa=theta
    hod=np.array([log_mcut,log_m1,sigma,alpha,kappa])
    nbins=fit_common['nmassbins'][0]

    if not(params['min'][0] < log_mcut < params['max'][0] and
           params['min'][1] < log_m1 < params['max'][1] and
           params['min'][2] < sigma < params['max'][2] and
           params['min'][3] < alpha < params['max'][3] and
           params['min'][4] < kappa < params['max'][4]):
        return -np.inf
    log_like=0.0
    N=hodmodel.num_gal(hod,massbin)
    wp2_obs=c2c.wp2_weighted(c2c.DD(hod,HH,HP,PP,wp2,fit_common,massbin),ar.RR(N,wp2),wp2)
    obs=wp2_obs
    for i in range(npoints):
        for j in range(npoints):
            log_like += (obs[i]-data[i])*invcov[i][j]*(obs[j]-data[j])
    return -0.5*log_like

def MCMC_hod(steps):
    BEfilename=chainsOut
    backend=emcee.backends.HDFBackend(BEfilename)
    backend.reset(16,5)

    with Pool() as pool: #parallel mcmc
        initial = params['init']
        ndim = len(initial)
        nwalkers = 16
        p0 = [np.array(initial)+0.1*np.random.randn(ndim) for i in range(nwalkers)]
        sampler = emcee.EnsembleSampler(nwalkers,ndim,log_probability,pool=pool,backend=backend)
        state=sampler.run_mcmc(p0,steps,progress=True)



if __name__=='__main__':
    
    if(dofitting):
        cov=np.loadtxt(covIn)
        invcov=np.linalg.inv(cov)
        if(wp2['usewp2']):
            data2=np.loadtxt(data2In)
            if(wp2['doPaircounts']):
                halocat=hodmodel.hal_cat(0,cosmo,pc_common,simulationIn,123)
                haloxyzn=np.float32(
                                np.stack(
                                (np.concatenate(halocat[0],axis=None),
                                np.concatenate(halocat[1],axis=None),
                                np.concatenate(halocat[2],axis=None),
                                np.concatenate(halocat[6],axis=None))))
                partxyzn=hodmodel.get_particle_catalog(halocat)
                massbin=pc.massbin_config(haloxyzn,partxyzn,fit_common,cosmo)
                pc_tmp=pc.rppipc(massbin,wp2,fit_common,pc_common,haloxyzn,partxyzn,paircountsOut)
                HH=pc_tmp[0]
                HP=pc_tmp[1]
                PP=pc_tmp[2]
            if(wp2['inputPaircounts']):
                HH=c2c.read_XX(paircountsIn,'HH',wp2['nrpbins'][0]*wp2['npibins'][0],fit_common['nmassbins'][0])
                HP=c2c.read_XY(paircountsIn,'HP',wp2['nrpbins'][0]*wp2['npibins'][0],fit_common['nmassbins'][0])
                PP=c2c.read_XX(paircountsIn,'PP',wp2['nrpbins'][0]*wp2['npibins'][0],fit_common['nmassbins'][0])
                massbin_mid=np.loadtxt(massbinmidIn)
                num_massbin_halo=np.loadtxt(nummassbinhaloIn)
                num_massbin_part=np.loadtxt(nummassbinpartIn)
                massbin=np.float32(np.stack((num_massbin_halo,num_massbin_part,massbin_mid)))
        if(xi3['usexi3']):
            data3=np.loadtxt(data3In)
            if(xi3['inputTrianglecounts']):
                HHH=c2c.read_XXX(triangleIn,'HHH',xi3['nsbins'][0],fit_common['nmassbins'][0])
                HPP=c2c.read_XYY(triangleIn,'HPP',xi3['nsbins'][0],fit_common['nmassbins'][0])
                PHH=c2c.read_XYY(triangleIn,'PHH',xi3['nsbins'][0],fit_common['nmassbins'][0])
                PPP=c2c.read_XXX(triangleIn,'PPP',xi3['nsbins'][0],fit_common['nmassbins'][0])
        npoints=wp2['nrpbins'][0]
        data=data2
        MCMC_hod(20000)
          

