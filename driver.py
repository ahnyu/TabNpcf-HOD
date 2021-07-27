###dependency###
import numpy as np
import pandas as pd
from multiprocessing import Pool
import time
###local file###
import ioconfig
from ioconfig import setParams,simParams,LRG,ELG,QSO,LXE,LXQ,EXQ,pathIn,binParams
import hodmodel
import emcee
import paircounts as pc
import counts2clustering as c2c
import analyticalRandom as ar
import loglike 


def MCMC_hod(steps):
    BEfilename='test.h5'
    backend=emcee.backends.HDFBackend(BEfilename)
    backend.reset(nparams_all*2,nparams_all)

    with Pool() as pool: #parallel mcmc
        initial = init_all
        ndim = len(initial)
        nwalkers = nparams_all*2
        p0 = [initial+width_all*np.random.randn(ndim) for i in range(nwalkers)]
        sampler = emcee.EnsembleSampler(nwalkers,ndim,loglike.loglike,
                args=(pairCounts,triCounts,setParams,tabParams,binParams,simParams,LRG,ELG,QSO,LXE,LXQ,EXQ,label_l,label_e,label_q,min_all,max_all,obs,obscov)
                ,pool=pool,backend=backend)
        state=sampler.run_mcmc(p0,steps,progress=True)



if __name__=='__main__':
    
    nmassbins=binParams['tab']['nmassbins']
    tabParams=np.empty((nmassbins,3))
    obs=np.empty((0))
    if(setParams['useWp']):
        wp_obs=np.loadtxt(pathIn['wpdata'],usecols=(1))
        rppibins=np.empty((2),dtype=int)
        obs=np.append(obs,wp_obs)
    if(setParams['useXil']):
        xil_obs=np.loadtxt(pathIn['xildata'])
        smubins=np.empty((2),dtype=int)
        obs=np.append(obs,xil_obs)
    if(setParams['useWp3']):
        wp3_obs=np.loadtxt(pathIn['wp3data'])
        obs=np.append(obs,wp3_obs)
    if(setParams['useXi3']):
        xi3_obs=np.loadtxt(pathIn['xi3data'])
        obs=np.append(obs,xi3_obs)
    obscov=np.loadtxt(pathIn['cov'])
    

    tabParams[:,0]=np.loadtxt(pathIn['tab']+'massBinMean.dat')
    tabParams[:,1]=np.loadtxt(pathIn['tab']+'numHaloBin.dat')
    tabParams[:,2]=np.loadtxt(pathIn['tab']+'numPartBin.dat')

    pairCounts={}
    triCounts={}
    if(setParams['useWp']):
        rppibins=np.array([binParams['rppi']['nrpbins'],binParams['rppi']['npibins']])
        npoints=rppibins[0]*rppibins[1]
        pairCounts['HH_rppi']=c2c.read_XX(pathIn['rppi'],'HH',npoints,nmassbins)
        pairCounts['HP_rppi']=c2c.read_XY(pathIn['rppi'],'HP',npoints,nmassbins)
        pairCounts['PP_rppi']=c2c.read_XX(pathIn['rppi'],'PP',npoints,nmassbins)
        print('read rppi HH HP PP finished')

    if(setParams['useXil']):
        smubins=np.array([binParams['smu']['nsbins'],binParams['smu']['nmubins']])
        npoints=smubins[0]*smubins[1]
        pairCounts['HH_smu']=c2c.read_XX(pathIn['smu'],'HH',npoints,nmassbins)
        pairCounts['HH_smu']=c2c.read_XY(pathIn['smu'],'HP',npoints,nmassbins)
        pairCounts['HH_smu']=c2c.read_XX(pathIn['smu'],'PP',npoints,nmassbins)
        print('read smu HH HP PP finished')

    if(setParams['useWp3']):
        triXYbins=binParams['triXY']['nsbins']
        triCounts['HHH_triXY']=c2c.read_XXX(pathIn['triXY'],'HHH',triXYbins,20)
        triCounts['HPP_triXY']=c2c.read_XYY(pathIn['triXY'],'HPP',triXYbins,20)
        triCounts['PHH_triXY']=c2c.read_XYY(pathIn['triXY'],'PHH',triXYbins,20)
        triCounts['PPP_triXY']=c2c.read_XXX(pathIn['triXY'],'PPP',triXYbins,20)
        print('read triXY HHH HPP PHH PPP finished')
        
    if(setParams['useXi3']):
        tri3Dbins=binParams['triXY']['nsbins']
        triCounts['HHH_tri']=c2c.read_XXX(pathIn['tri3D'],'HHH',tri3Dbins,20)
        triCounts['HPP_tri']=c2c.read_XYY(pathIn['tri3D'],'HPP',tri3Dbins,20)
        triCounts['PHH_tri']=c2c.read_XYY(pathIn['tri3D'],'PHH',tri3Dbins,20)
        triCounts['PPP_tri']=c2c.read_XXX(pathIn['tri3D'],'HHH',tri3Dbins,20)
        print('read tri3D HHH HPP PHH PPP finished')

    npoints_all=0
    n2pt_all=0
    n2pt_l=0;n2pt_e=0;n2pt_q=0;
    n2pt_lxe=0;n2pt_lxq=0;n2pt_exq=0;
    n3pt_all=0
    n3pt_l=0;n3pt_e=0;n3pt_q=0;
    nparams_all=0
    init_all=np.empty((0))
    min_all=np.empty((0))
    max_all=np.empty((0))
    width_all=np.empty((0))
    label_l=np.empty((0),dtype=str)
    label_e=np.empty((0),dtype=str)
    label_q=np.empty((0),dtype=str)

    if(setParams['useLRG']):
        if(setParams['useWp']):
            n2pt_l=LRG['wpidx']['max']-LRG['wpidx']['min']
            print('using LRG wp')
        elif(setParams['useXil']):
            n2pt_l=LRG['xi0idx']['max']-LRG['xi0idx']['min']+LRG['xi2idx']['max']-LRG['xi2idx']['min']
            print('using LRG xil')
        if(setParams['useWp3']):
            n3pt_l=LRG['wp3idx']['max']-LRG['wp3idx']['min']
            print('using LRG wp3')
        elif(setParams['useXi3']):
            n3pt_l=LRG['xi3idx']['max']-LRG['xi3idx']['min']
            print('using LRG xi3')
        n2pt_all+=n2pt_l
        n3pt_all+=n3pt_l

        if(LRG['model']['cent']==1):
            if(setParams['useELG'] or setParams['useQSO']):
                label_l=np.append(label_l,['logMcut','sigma','pmax'])
            else:
                label_l=np.append(label_l,['logMcut','sigma'])
        if(LRG['model']['sate']==1):
            label_l=np.append(label_l,['logM1','alpha','kappa'])
        nparams_l=len(label_l)

        df_l=pd.DataFrame.from_dict(LRG['params'])
        init_all=np.append(init_all,df_l.loc['init',label_l].to_numpy()) 
        min_all=np.append(min_all,df_l.loc['min',label_l].to_numpy()) 
        max_all=np.append(max_all,df_l.loc['max',label_l].to_numpy()) 
        width_all=np.append(width_all,df_l.loc['width',label_l].to_numpy()) 

        nparams_all+=nparams_l
        print('num of 2pcf points = ',n2pt_l,', num of 3pcf points = ',n3pt_l)
        print('LRG HOD model, num of params =',nparams_l)

    if(setParams['useELG']):
        if(setParams['useWp']):
            n2pt_e=ELG['wpidx']['max']-ELG['wpidx']['min']
            print('using ELG wp')
        elif(setParams['useXil']):
            n2pt_e=ELG['xi0idx']['max']-ELG['xi0idx']['min']+ELG['xi2idx']['max']-ELG['xi2idx']['min']
            print('using ELG xil')
        if(setParams['useWp3']):
            n3pt_e=ELG['wp3idx']['max']-ELG['wp3idx']['min']
            print('using ELG wp3')
        elif(setParams['useXi3']):
            n3pt_e=ELG['xi3idx']['max']-ELG['xi3idx']['min']
            print('using ELG xi3')
        n2pt_all+=n2pt_e
        n3pt_all+=n3pt_e


        if(ELG['model']['cent']==1):
            label_e=np.append(label_e,['logMcut','sigma','pmax'])
            if(ELG['model']['sate']==1):
                label_e=np.append(label_e,['logM1','alpha','kappa'])
        elif(ELG['model']['cent']==2):
            label_e=np.append(label_e,['logMcut','sigma','pmax','gamma','Q'])
            if(ELG['model']['sate']==1):
                label_e=np.append(label_e,['logM1','alpha'])
        nparams_e=len(label_e)

        df_e=pd.DataFrame.from_dict(ELG['params'])
        init_all=np.append(init_all,df_e.loc['init',label_e].to_numpy()) 
        min_all=np.append(min_all,df_e.loc['min',label_e].to_numpy()) 
        max_all=np.append(max_all,df_e.loc['max',label_e].to_numpy()) 
        width_all=np.append(width_all,df_e.loc['width',label_e].to_numpy()) 

        nparams_all+=nparams_e
        print('num of 2pcf points = ',n2pt_e,', num of 3pcf points = ',n3pt_e)
        print('ELG HOD model, num of params =',nparams_e)


    if(setParams['useQSO']):
        if(setParams['useWp']):
            n2pt_q=QSO['wpidx']['max']-QSO['wpidx']['min']
            print('using QSO wp')
        elif(setParams['useXil']):
            n2pt_q=QSO['xi0idx']['max']-QSO['xi0idx']['min']+QSO['xi2idx']['max']-QSO['xi2idx']['min']
            print('using QSO xil')
        if(setParams['useWp3']):
            n3pt_q=QSO['wp3idx']['max']-QSO['wp3idx']['min']
            print('using QSO wp3')
        elif(setParams['useXi3']):
            n3pt_q=QSO['xi3idx']['max']-QSO['xi3idx']['min']
            print('using QSO xi3')
        n2pt_all+=n2pt_q
        n3pt_all+=n3pt_q


        if(QSO['model']['cent']==1):
            label_q=np.append(label_q,['logMcut','sigma','pmax'])
        if(QSO['model']['sate']==1):
            label_q=np.append(label_q,['logM1','alpha','kappa'])
        elif(QSO['model']['sate']==2):
            label_q=np.append(label_q,['logM1','alpha','logMmin'])

        nparams_q=len(label_q)

        df_q=pd.DataFrame.from_dict(QSO['params'])
        init_all=np.append(init_all,df_q.loc['init',label_q].to_numpy()) 
        min_all=np.append(min_all,df_q.loc['min',label_q].to_numpy()) 
        max_all=np.append(max_all,df_q.loc['max',label_q].to_numpy()) 
        width_all=np.append(width_all,df_q.loc['width',label_q].to_numpy()) 

        nparams_all+=nparams_q
        print('num of 2pcf points = ',n2pt_q,', num of 3pcf points = ',n3pt_q)
        print('QSO HOD model, num of params =',nparams_q)

    if(setParams['useLXE']):
        n2pt_lxe=LXE['wpidx']['max']-LXE['wpidx']['min']
        n2pt_all+=n2pt_lxe
        print('using LXE wp, num of points =',n2pt_lxe)
    if(setParams['useLXQ']):
        n2pt_lxq=LXQ['wpidx']['max']-LXQ['wpidx']['min']
        n2pt_all+=n2pt_lxq
        print('using LXQ wp, num of points =',n2pt_lxq)
    if(setParams['useEXQ']):
        n2pt_exq=EXQ['wpidx']['max']-EXQ['wpidx']['min']
        n2pt_all+=n2pt_exq
        print('using EXQ wp, num of points =',n2pt_exq)

    npoints_all=n2pt_all+n3pt_all

    if(npoints_all==len(obs)):
        print('total points =',npoints_all,', 2pcf points = ',n2pt_all,', 3pcf points = ',n3pt_all)
    else:
        print('points from theory and data do not match please check')



    st=time.time()
#    print(log_probability(init_all))



    print(loglike.loglike(init_all,pairCounts,triCounts,setParams,tabParams,binParams,simParams,LRG,ELG,QSO,LXE,LXQ,EXQ,label_l,label_e,label_q,min_all,max_all,obs,obscov))
        
    end=time.time()


    print(end-st)

    MCMC_hod(20000)
          

