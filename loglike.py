import numpy as np
import math
import counts2clustering as c2c
import analyticalRandom as ar
import hodmodel

def loglike(theta,pairCounts,triCounts,setParams,tabParams,binParams,simParams,LRG,ELG,QSO,LXE,LXQ,EXQ,label_l,label_e,label_q,min_all,max_all,obs,obscov):

    total_pmax=0.
    idx=0
    if(setParams['useLRG']):
        params_l={}
        for key in label_l:
            params_l[key]=theta[idx]
            idx+=1
        if not(setParams['useELG'] or setParams['useLRG']):
            params_l['pmax']=1.0
        total_pmax+=params_l['pmax']
    if(setParams['useELG']):
        params_e={}
        for key in label_e:
            params_e[key]=theta[idx]
            idx+=1
        if(ELG['model']['cent']==2):
            params_e['kappa']=1.0
            params_e['maxpdf']=hodmodel.findmaxpdf(params_e['logMcut'],params_e['sigma'],params_e['gamma'])
        total_pmax+=params_e['pmax']
    if(setParams['useQSO']):
        params_q={}
        for key in label_q:
            params_q[key]=theta[idx]
            idx+=1
        total_pmax+=params_q['pmax']


    if(setParams['useWp']):
        HH_rppi=pairCounts['HH_rppi']
        HP_rppi=pairCounts['HP_rppi']
        PP_rppi=pairCounts['PP_rppi']
        rppibins=np.array([binParams['rppi']['nrpbins'],binParams['rppi']['npibins']],dtype=int)
    if(setParams['useXil']):
        HH_smu=pairCounts['HH_smu']
        HP_smu=pairCounts['HP_smu']
        PP_smu=pairCounts['PP_smu']
        smubins=np.array([binParams['smu']['nsbins'],binParams['smu']['nmubins']],dtype=int)
    if(setParams['useWp3']):
        HHH_triXY=triCounts['HHH_triXY']
        HPP_triXY=triCounts['HPP_triXY']
        PHH_triXY=triCounts['PHH_triXY']
        PPP_triXY=triCounts['PPP_triXY']
    if(setParams['useXi3']):
        HHH_tri=triCounts['HHH_tri3D']
        HPP_tri=triCounts['HPP_tri3D']
        PHH_tri=triCounts['PHH_tri3D']
        PPP_tri=triCounts['PPP_tri3D']



    if not(np.all(theta<max_all) and np.all(theta>min_all)):
        return -np.inf
    if (total_pmax>1):
        return -np.inf

    wp_the=np.empty((0))
    xil_the=np.empty((0))
    wp3_the=np.empty((0))
    xi3_the=np.empty((0))

    if(setParams['useLRG']):
        N_lrg=hodmodel.num_lrg(params_l,tabParams[:,0],tabParams[:,1],LRG['model'])
        if(setParams['useWp']):
            wp_lrg=c2c.wp2(c2c.DD_lrg(params_l,HH_rppi,HP_rppi,PP_rppi,rppibins,tabParams,LRG['model'])\
                    ,ar.RR(N_lrg,N_lrg,binParams,simParams['boxsize']),rppibins)\
                    [LRG['wpidx']['min']:LRG['wpidx']['max']]
            wp_the=np.append(wp_the,wp_lrg)
    if(setParams['useELG']):
        N_elg=hodmodel.num_elg(params_e,tabParams[:,0],tabParams[:,1],ELG['model'])
        if(setParams['useWp']):
            wp_elg=c2c.wp2(c2c.DD_elg(params_e,HH_rppi,HP_rppi,PP_rppi,rppibins,tabParams,ELG['model'])
                    ,ar.RR(N_elg,N_elg,binParams,simParams['boxsize']),rppibins)\
            [ELG['wpidx']['min']:ELG['wpidx']['max']]
            wp_the=np.append(wp_the,wp_elg)
    if(setParams['useQSO']):
        N_qso=hodmodel.num_qso(params_q,tabParams[:,0],tabParams[:,1],QSO['model'])
        if(setParams['useWp']):
            wp_qso=c2c.wp2(c2c.DD_qso(params_q,HH_rppi,HP_rppi,PP_rppi,rppibins,tabParams,QSO['model'])
                    ,ar.RR(N_qso,N_qso,binParams,simParams['boxsize']),rppibins)\
            [QSO['wpidx']['min']:QSO['wpidx']['max']]
            wp_the=np.append(wp_the,wp_qso)

    if(setParams['useLXE']):
        if(setParams['useWp']):
            wp_lxe=c2c.wp2(c2c.DD_lxe(params_l,params_e,HH_rppi,HP_rppi,PP_rppi,rppibins,tabParams,LRG['model'],ELG['model'])
                    ,ar.RR(N_lrg,N_elg,binParams,simParams['boxsize']),rppibins)\
            [LXE['wpidx']['min']:LXE['wpidx']['max']]
            wp_the=np.append(wp_the,wp_lxe)
    if(setParams['useLXQ']):
        if(setParams['useWp']):
            wp_lxq=c2c.wp2(c2c.DD_lxq(params_l,params_q,HH_rppi,HP_rppi,PP_rppi,rppibins,tabParams,LRG['model'],QSO['model'])
                    ,ar.RR(N_lrg,N_qso,binParams,simParams['boxsize']),rppibins)\
            [LXQ['wpidx']['min']:LXQ['wpidx']['max']]
            wp_the=np.append(wp_the,wp_lxq)
    if(setParams['useEXQ']):
        if(setParams['useWp']):
            wp_exq=c2c.wp2(c2c.DD_exq(params_e,params_q,HH_rppi,HP_rppi,PP_rppi,rppibins,tabParams,ELG['model'],QSO['model'])
                    ,ar.RR(N_elg,N_qso,binParams,simParams['boxsize']),rppibins)\
            [EXQ['wpidx']['min']:EXQ['wpidx']['max']]
            wp_the=np.append(wp_the,wp_exq)

    the=np.hstack((wp_the,xil_the,wp3_the,xi3_the))

    diff=the-obs
#    print(the)
    log_like=-0.5*np.dot(diff,np.linalg.solve(obscov,diff))

    if(math.isnan(log_like)):
        return -np.inf


    return log_like
