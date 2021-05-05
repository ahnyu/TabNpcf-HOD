import argparse
from configobj import ConfigObj
import numpy as np

parser=argparse.ArgumentParser()
parser.add_argument('--inifile','-ini',help='set config file path')
args=parser.parse_args()



config=ConfigObj(args.inifile)

settings=config['settings']
dofitting=settings.as_bool('doFitting')
HODcatalog=settings.as_bool('HODcatalog')


if(dofitting):
    fit_dic=config['fit']
    ###params###
    params_dic=fit_dic['params']
    params_names=params_dic.keys()
    params=np.empty(len(params_names),dtype={'names':('min','max','init'),'formats':('float32','float32','float32')})
    for i in range(len(params_names)):
        params['min'][i]=params_dic[params_names[i]]['min']
        params['max'][i]=params_dic[params_names[i]]['max']
        params['init'][i]=params_dic[params_names[i]]['init']

    ###use2pt###
    use2pt_dic=fit_dic['use2pt']
    wp2_dic=use2pt_dic['wp2']
    wp2_names=wp2_dic.keys()
    wp2=np.zeros(1,dtype={'names':wp2_names,'formats':[np.bool_,np.bool_,np.bool_,'f4','f4',np.int,np.int,'f4']})
    wp2['usewp2']=wp2_dic.as_bool('usewp2')
    if(wp2['usewp2']):
        for i in range(1,len(wp2_names)):
            if(wp2[wp2_names[i]].dtype==bool):
                wp2[wp2_names[i]]=wp2_dic.as_bool(wp2_names[i])
            else:
                wp2[wp2_names[i]]=wp2_dic[wp2_names[i]]
            print(wp2_names[i],wp2[wp2_names[i]])
        ###paircounts###
        if(wp2['doPaircounts']):
            pc_dic=config['paircounts']
            ###paircounts common###
            pc_common_dic=pc_dic['common']
            pc_common_names=pc_common_dic.keys()
            pc_common=np.empty(1,dtype={'names':pc_common_names,'formats':[np.int,'f4',np.bool_]})
            for i in range(len(pc_common_names)):
                pc_common[pc_common_names[i]]=pc_common_dic[pc_common_names[i]]
            ###paircountsOut###
            if(pc_common['savepc']):
                paircountsOut=pc_dic['paircountsOut']
            ###simulationIn###
            simulationIn=pc_dic['simulationIn']
        if(wp2['inputPaircounts']):
            paircountsIn=config['readCounts']['paircountsIn']
    ###use3pt###
    use3pt_dic=fit_dic['use3pt']
    xi3_dic=use3pt_dic['xi3']
    xi3_names=xi3_dic.keys()
    xi3=np.zeros(1,dtype={'names':xi3_names,'formats':[np.bool_,np.bool_,np.int]})
    xi3['usexi3']=xi3_dic.as_bool('usexi3')
    if(xi3['usexi3']):
        for i in range(1,len(xi3_names)):
            if(xi3[xi3_names[i]].dtype==bool):
                xi3[xi3_names[i]]=xi3_dic.as_bool(xi3_names[i])
            else:
                xi3[xi3_names[i]]=xi3_dic[xi3_names[i]]
        if(xi3['inputTrianglecounts']):
            triangleIn=config['readCounts']['trianglecountsIn']
    ###fit common###
    fit_common_dic=fit_dic['common']
    fit_common_names=fit_common_dic.keys()
    fit_common=np.ones(1,dtype={'names':fit_common_names,'formats':[np.int,np.int]})
    for i in range(len(fit_common_names)):
        if(fit_common[fit_common_names[i]].dtype==bool):
            fit_common[fit_common_names[i]]=fit_common_dic.as_bool(fit_common_names[i])
        else:
            fit_common[fit_common_names[i]]=fit_common_dic[fit_common_names[i]]
    ###chainsOut###
    chainsOut=fit_dic['chainsOut']
    ###covIn###
    covIn=fit_dic['covIn']
    ###data2In###
    data2In=fit_dic['data2In']
    ###data3In###
    data3In=fit_dic['data3In']
###massbin setting###
if(wp2['inputPaircounts'] or xi3['inputTrianglecounts']):
    massbinmidIn=config['readCounts']['massbinmidIn']
    nummassbinhaloIn=config['readCounts']['nummassbinhaloIn']
    nummassbinpartIn=config['readCounts']['nummassbinpartIn']
###HODcat###
if(HODcatalog):
    HODcat_dic=config['HODcat']
    HODcat_names=HODcat_dic.keys()
    HODcat=np.empty(1,dtype={'names':HODcat_names,'formats':['f4','f4','f4','f4','f4']})
    for i in range(len(HODcat_names)):
        HODcat[HODcat_names[i]]=HODcat_dic[HODcat_names[i]]
    
###cosmology###
cosmo_dic=config['cosmology']
cosmo_names=cosmo_dic.keys()
cosmo=np.empty(1,dtype={'names':cosmo_names,'formats':['f4','f4','f8']})
for i in range(len(cosmo_names)):
    cosmo[cosmo_names[i]]=cosmo_dic[cosmo_names[i]]

