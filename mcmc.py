import numpy as np
from multiprocessing import Pool
import emcee
import hodmodel 

def log_probability(theta,params,data,invcov,fit_common,npoints,HH,HP,PP):
    log_mcut,log_m1,sigma,alpha,kappa=theta
    hod=np.array([log_mcut,log_m1,sigma,alpha,kappa])
    nbins=fit_common['nmassbins']

    if not(params[0]['min'] < log_mcut < prior[0]['max'] and
           params[1]['min'] < log_m1 < prior[1]['max'] and
           params[2]['min'] < sigma < prior[2]['max'] and
           params[3]['min'] < alpha < prior[3]['max'] and
           params[4]['min'] < kappa < prior[4]['max']):
        return -np.inf
    log_like=0.0
    N=num_gal(hod)
    wp2_obs=wp2(DD(hod,nbins),RR(N))
    obs=wp2_obs


    for i in range(npoints):
        for j in range(npoints):
            log_like += (obs[i]-data[i])*invcov[i][j]*(obs[j]-data[j])
    return -0.5*log_like


###MCMC
def MCMC_hod(chainsOut,params,data,invcov,npoints)
    BEfilename=chainsOut
    backend=emcee.backends.HDFBackend(BEfilename)
    backend.reset(16,5)

    with Pool() as pool: #parallel mcmc
        initial = params['init']
        ndim = len(initial)
        nwalkers = 16
        p0 = [np.array(initial)+0.1*np.random.randn(ndim) for i in range(nwalkers)]
        sampler = emcee.EnsembleSampler(nwalkers,ndim,log_probability,args(params,data,invcov,fit_common,npoints),pool=pool,backend=backend)
        state=sampler.run_mcmc(p0,steps,progress=True)
