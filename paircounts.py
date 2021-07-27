import numpy as np
from configobj import ConfigObj
from Corrfunc.theory.DDrppi import DDrppi
from numba import jit,njit
from astropy.table import Table
from argparse import ArgumentParser

def read_downsamp(root,nfiles):
    haloscat=[]
    partscat=[]
    halos=Table()
    parts=Table()
    label=['pos','Mhalo']
    for i in range(nfiles):
        haloscat.append(Table.read(root+'halos_%02d.hdf5' %i,path='data'))
        partscat.append(Table.read(root+'parts_%02d.hdf5' %i,path='data'))
    htmppos=np.empty((0,3),dtype=np.float32)
    ptmppos=np.empty((0,3),dtype=np.float32)
    
    htmpM=np.empty((0))
    ptmpM=np.empty((0))
    for j in range(nfiles):
        htmppos=np.append(htmppos,haloscat[j][label[0]],axis=0)
        ptmppos=np.append(ptmppos,partscat[j][label[0]],axis=0)
        htmpM=np.append(htmpM,haloscat[j][label[1]])
        ptmpM=np.append(ptmpM,partscat[j][label[1]])
        
    halos[label[0]]=htmppos
    halos[label[1]]=htmpM
    
    parts[label[0]]=ptmppos
    parts[label[1]]=ptmpM

    return halos,parts

@njit(parallel=True)
def tabulation(hmass,phmass,Mpart,nbins):
    massbin_edge=np.linspace(np.log10(np.min(hmass)),np.log10(np.max(hmass)+Mpart),nbins+1)
    massbin=np.zeros(nbins)
    num_halo=np.full(nbins,0)
    num_part=np.full(nbins,0)
    mask_halo=np.full(len(hmass),-1)
    mask_part=np.full(len(phmass),-1)
    for i in range(nbins):
        submask_halo=(np.log10(hmass)>=massbin_edge[i])&(np.log10(hmass)<massbin_edge[i+1])
        submask_part=(np.log10(phmass)>=massbin_edge[i])&(np.log10(phmass)<massbin_edge[i+1])
        mask_halo[submask_halo]=i
        mask_part[submask_part]=i
        num_halo[i]=len(mask_halo[submask_halo])
        num_part[i]=len(mask_part[submask_part])
        massbin[i]=np.log10(np.mean(hmass[submask_halo]))
    return num_halo,num_part,massbin_edge,massbin,mask_halo,mask_part

def rppipc(massmask_halo,massmask_part,halo,part,nthreads,pimax,rpbins,boxsize,nbins,outroot):

    halox=halo['pos'][:,0]
    haloy=halo['pos'][:,1]
    haloz=halo['pos'][:,2]

    partx=part['pos'][:,0]
    party=part['pos'][:,1]
    partz=part['pos'][:,2]


    print('start tabulation paircounting (rppi)')
    for i in range(nbins):
        for j in range(nbins):
            print('working on massbin %s' %i+' and %s' %j)
            if(i==j):
                autocorr=1
                HH=DDrppi(autocorr,nthreads,pimax,rpbins,
                        halox[massmask_halo==i],haloy[massmask_halo==i],
                        haloz[massmask_halo==i],boxsize=boxsize)
                PP=DDrppi(autocorr,nthreads,pimax,rpbins,
                        partx[massmask_part==i],party[massmask_part==i],
                        partz[massmask_part==i],boxsize=boxsize)

                np.savetxt(outroot+'HH_%02d' %i +'_%02d.dat' %j, HH)
                np.savetxt(outroot+'PP_%02d' %i +'_%02d.dat' %j, PP)
            elif(i>j):
                autocorr=0
                HH=DDrppi(autocorr,nthreads,pimax,rpbins,
                        X1=halox[massmask_halo==i],Y1=haloy[massmask_halo==i],
                        Z1=haloz[massmask_halo==i],X2=halox[massmask_halo==j],
                        Y2=haloy[massmask_halo==j],Z2=haloz[massmask_halo==j],
                        boxsize=boxsize)
                PP=DDrppi(autocorr,nthreads,pimax,rpbins,
                        X1=partx[massmask_part==i],Y1=party[massmask_part==i],
                        Z1=partz[massmask_part==i],X2=partx[massmask_part==j],
                        Y2=party[massmask_part==j],Z2=partz[massmask_part==j],
                        boxsize=boxsize)

                np.savetxt(outroot+'HH_%02d' %i +'_%02d.dat' %j, HH)
                np.savetxt(outroot+'PP_%02d' %i +'_%02d.dat' %j, PP)
            autocorr=0
            HP=DDrppi(autocorr,nthreads,pimax,rpbins,
                        X1=halox[massmask_halo==i],Y1=haloy[massmask_halo==i],
                        Z1=haloz[massmask_halo==i],X2=partx[massmask_part==j],
                        Y2=party[massmask_part==j],Z2=partz[massmask_part==j],
                        boxsize=boxsize)

            np.savetxt(outroot+'HP_%02d' %i +'_%02d.dat' %j, HP)
    return HH,HP,PP


def main():
    parser = ArgumentParser(description='tabulation paricount')
    parser.add_argument('--inifile','-ini',help ='config file for downsample')
    args=parser.parse_args()

    config=ConfigObj(args.inifile)

    clustering=config['clustering']
    downcatpath=config['downcatpath']
    outpath=config['outpath']
    nfiles=config.as_int('nfiles')
    nrpbins=config.as_int('nrpbins')
    rpmin=config.as_float('rpmin')
    rpmax=config.as_float('rpmax')
    pimax=config.as_int('pimax')
    nthreads=config.as_int('nthreads')
    boxsize=config.as_float('boxsize')
    nmassbins=config.as_int('nmassbins')
    Mpart=config.as_float('particlemass')

    rpbins=np.logspace(rpmin,rpmax,nrpbins+1)

    print('reading downsample')
    halosall,partsall=read_downsamp(downcatpath,nfiles)
    print('reading done, tabulating')
    num_halo,num_part,bin_edge,bin_mid,mask_halo,mask_part=tabulation(halosall['Mhalo'],partsall['Mhalo'],Mpart,nmassbins)
    np.savetxt(outpath+'numHaloBin.dat',num_halo,fmt='%i')
    np.savetxt(outpath+'numPartBin.dat',num_part,fmt='%i')
    np.savetxt(outpath+'massBinEdge.dat',bin_edge)
    np.savetxt(outpath+'massBinMean.dat',bin_mid)
    print('tabulating done, paircounting')
    if(clustering=='rppi'):
        rppipc(mask_halo,mask_part,halosall,partsall,nthreads,pimax,rpbins,boxsize,nmassbins,outpath)
    else:
        print('other clustering type to be update')
    return 0

if __name__ == '__main__':
    main()

