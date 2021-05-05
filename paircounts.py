import numpy as np
from Corrfunc.theory.DDrppi import DDrppi

def massbin_config(halo,part,fit_common,cosmo):
    nbins=fit_common['nmassbins'][0]
    mass_particle=cosmo['particleMass']
    halon=halo[3]
    partn=part[3]
    num_massbin_halo=[]
    massmask_halo=[]
    num_massbin_part=[]
    massmask_part=[]
    massbin_mid=[]
    massbinlog=np.linspace(np.log10(np.min(halon)*mass_particle),np.log10((np.max(halon)+1)*mass_particle),nbins+1)
    for i in range(nbins):
        mask_halo=(np.log10(halon*mass_particle)>=massbinlog[i])&(np.log10(halon*mass_particle)<massbinlog[i+1])
        mask_part=(np.log10(partn*mass_particle)>=massbinlog[i])&(np.log10(partn*mass_particle)<massbinlog[i+1])
        massmask_halo.append(mask_halo)
        massmask_part.append(mask_part)
        num_massbin_halo.append(len(halon[mask_halo]))
        num_massbin_part.append(len(partn[mask_part]))    
        massbin_mid.append(np.mean(halon[mask_halo]))
    num_massbin_halo=np.array(num_massbin_halo)
    num_massbin_part=np.array(num_massbin_part)
    massmask_halo=np.array(massmask_halo)
    massmask_part=np.array(massmask_part)
    massbin_mid=np.array(massbin_mid)
    
    return num_massbin_halo,num_massbin_part,massbin_mid,massmask_halo,massmask_part

def rppipc(massbin,wp2,fit_common,pc_common,halo,part,paircountsOut):
    
    nrpbins=wp2['nrpbins'][0]
    npibins=wp2['npibins'][0]
    rpbins=np.logspace(wp2['logrpmin'][0],wp2['logrpmax'][0],nrpbins+1)
    pimax=wp2['pimax'][0]
    nthreads=pc_common['nthreads'][0]
    boxsize=pc_common['boxsize'][0]
    savepc=pc_common['savepc'][0]
    nbins=fit_common['nmassbins'][0]
    massmask_halo=massbin[3]
    massmask_part=massbin[4]
    halox=halo[0]
    haloy=halo[1]
    haloz=halo[2]
    partx=part[0]
    party=part[1]
    partz=part[2]
    outroot=paircountsOut
    print(halox.dtype)
    print(partx.dtype)

    HHall=np.empty((nrpbins*npibins,nbins,nbins))
    HPall=np.empty((nrpbins*npibins,nbins,nbins))
    PPall=np.empty((nrpbins*npibins,nbins,nbins))
    for i in range(nbins):
        for j in range(nbins):
            print('working on i=%s' %i+' j=%s' %j)
            if(i==j):
                autocorr=1
                HH=DDrppi(autocorr,nthreads,pimax,rpbins,
                        halox[massmask_halo[i]],haloy[massmask_halo[i]],
                        haloz[massmask_halo[i]],boxsize=boxsize)
                PP=DDrppi(autocorr,nthreads,pimax,rpbins,
                        partx[massmask_part[i]],party[massmask_part[i]],
                        partz[massmask_part[i]],boxsize=boxsize)
                HHall[:,i,j]=HH['npairs']
                PPall[:,i,j]=PP['npairs']

                if(savepc):
                    np.savetxt(outroot+'HH_%02d' %i +'_%02d.dat' %j, HH)
                    np.savetxt(outroot+'PP_%02d' %i +'_%02d.dat' %j, PP)
            elif(i>j):
                autocorr=0
                HH=DDrppi(autocorr,nthreads,pimax,rpbins,
                        X1=halox[massmask_halo[i]],Y1=haloy[massmask_halo[i]],
                        Z1=haloz[massmask_halo[i]],X2=halox[massmask_halo[j]],
                        Y2=haloy[massmask_halo[j]],Z2=haloz[massmask_halo[j]],
                        boxsize=boxsize)
                PP=DDrppi(autocorr,nthreads,pimax,rpbins,
                        X1=partx[massmask_part[i]],Y1=party[massmask_part[i]],
                        Z1=partz[massmask_part[i]],X2=partx[massmask_part[j]],
                        Y2=party[massmask_part[j]],Z2=partz[massmask_part[j]],
                        boxsize=boxsize)
                HHall[:,i,j]=HH['npairs']
                HHall[:,j,i]=HH['npairs']
                PPall[:,i,j]=PP['npairs']
                PPall[:,j,i]=PP['npairs']

                if(savepc):
                    np.savetxt(outroot+'HH_%02d' %i +'_%02d.dat' %j, HH)
                    np.savetxt(outroot+'PP_%02d' %i +'_%02d.dat' %j, PP)
            autocorr=0
            HP=DDrppi(autocorr,nthreads,pimax,rpbins,
                        X1=halox[massmask_halo[i]],Y1=haloy[massmask_halo[i]],
                        Z1=haloz[massmask_halo[i]],X2=partx[massmask_part[j]],
                        Y2=party[massmask_part[j]],Z2=partz[massmask_part[j]],
                        boxsize=boxsize)
            HPall[:,i,j]=HP['npairs']
            if(savepc):
                np.savetxt(outroot+'HP_%02d' %i +'_%02d.dat' %j, HP)
    return HHall,HPall,PPall




