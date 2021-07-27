from abacusnbody.data.compaso_halo_catalog import CompaSOHaloCatalog
from astropy.table import Table
from argparse import ArgumentParser
from configobj import ConfigObj
import glob
import numpy as np
import gc

def downhalo(x):
    return 1.0/(1.0 + 10*np.exp(-(x - 11.2)*25))

def prep_halo_part(halos,parts,frac,rsd,step,seed,simParams):
    halosout=Table()
    partsout=Table()
    if(step):
        down_halo=downhalo(np.log10(halos['N']*simParams['Mpart']))
        cut_mask=(halos['N']>simParams['Ncut'])&(np.random.random(len(halos))<down_halo)
    else:
        cut_mask=(halos['N']>simParams['Ncut'])
    halos_cut=halos[cut_mask]
    halos_cut_pidx=halos['npstartA'][cut_mask]
    halos_cut_pnum=halos['npoutA'][cut_mask]
    halos_cut_pnum_new=np.full(len(halos_cut),0)
    parts_mask=np.full(len(parts),0)
    parts_halomass=np.full(len(parts),-1.)
    parts_ind=np.full(len(parts),-1)
    for i in np.arange(len(halos_cut)):
        if(i%5000000==0):
            print(i)
        tmp_mask=np.random.binomial(n=1,p=frac,size=halos_cut_pnum[i])
        halos_cut_pnum_new[i]=np.sum(tmp_mask)
        parts_mask[halos_cut_pidx[i]:halos_cut_pidx[i]+halos_cut_pnum[i]]=tmp_mask
        parts_halomass[halos_cut_pidx[i]:halos_cut_pidx[i]+halos_cut_pnum[i]]=halos_cut['N'][i]*simParams['Mpart']
        parts_ind[halos_cut_pidx[i]:halos_cut_pidx[i]+halos_cut_pnum[i]]=i

    halosout['Mhalo']=halos_cut['N']*simParams['Mpart']
    halosout['pos']=halos_cut['x_L2com']+1000.
    partsout['Mhalo']=parts_halomass[parts_mask==1]
    partsout['pos']=parts['pos'][parts_mask==1]+1000.

    if(rsd):
        halos_ztmp=halos_cut['x_L2com'][:,2]+1000.
        parts_ztmp=parts['pos'][:,2][parts_mask==1]+1000.
        halosout['pos'][:,2]=np.mod(halos_ztmp+halos_cut['v_L2com'][:,2]*simParams['rsd'],simParams['boxsize'])
        partsout['pos'][:,2]=np.mod(parts_ztmp+parts['vel'][:,2][parts_mask==1]*simParams['rsd'],simParams['boxsize'])


    np.random.seed(seed)
    halosout['randoms']=np.random.random(len(halosout))
    partsout['randoms']=np.random.random(len(partsout))
    halosout['nparts']=halos_cut_pnum_new
    partsout['hidx']=parts_ind[parts_mask==1]
    return halosout,partsout

def main():
    parser = ArgumentParser(description='prepare downsample abacus box for tabulation')
    parser.add_argument('--inifile','-ini',help ='config file for downsample')
    args=parser.parse_args()

    config=ConfigObj(args.inifile)

    abacusbox=config['abacusbox']
    redshift=config['redshift']
    outpath=config['outpath']
    frac=config.as_float('partfrac')
    seed=config.as_int('randseed')
    mcut=config.as_float('masscut')
    step=config.as_bool('halostep')
    rsd=config.as_bool('usersd')

    print(abacusbox,redshift,outpath,frac,seed,mcut,step,rsd)

    
    abacuscat=sorted(glob.glob('/global/cfs/cdirs/desi/cosmosim/Abacus/'+abacusbox+'/halos/'+redshift+'/halo_info/*.asdf'))

    testcat=CompaSOHaloCatalog(abacuscat[0],fields=['N'],cleaned=True)

    simParams={}
    simParams['redshift']=testcat.header['Redshift']
    simParams['h']=testcat.header['H0']
    simParams['boxsize']=testcat.header['BoxSize']
    simParams['rsd']=1/(testcat.header['VelZSpace_to_kms']/testcat.header['BoxSize'])
    simParams['Mpart']=testcat.header['ParticleMassHMsun']
    simParams['Ncut']=mcut/simParams['Mpart']
    print(simParams)
    print('Ncut =',simParams['Ncut'])

    for islab in range(9):   
    	print('working on slab',islab)
    	
    	if(islab==8):
    	    cat=CompaSOHaloCatalog(abacuscat[islab*4:], fields=['N','x_L2com','v_L2com','npstartA','npoutA'],subsamples=dict(A=True,rv=True),cleaned=True)
    	else:
    	    cat=CompaSOHaloCatalog(abacuscat[islab*4:4+islab*4], fields=['N','x_L2com','v_L2com','npstartA','npoutA'],subsamples=dict(A=True,rv=True),cleaned=True)
    	    
    	print('finish reading slab')
    	haloslab,partslab=prep_halo_part(cat.halos,cat.subsamples,frac,rsd,step,seed+islab,simParams)
    	print('finish prep')
    	haloslab.write(outpath+'halos_%02d.hdf5' %islab,path='data')
    	partslab.write(outpath+'parts_%02d.hdf5' %islab,path='data')
    	del cat
    	del haloslab
    	del partslab
    	gc.collect()

    return 0

if __name__ == '__main__':
    main()



