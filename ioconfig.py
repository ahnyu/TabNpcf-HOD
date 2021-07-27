import argparse
from configobj import ConfigObj
import numpy as np

parser=argparse.ArgumentParser()
parser.add_argument('--inifile','-ini',help='set config file path')
args=parser.parse_args()



config=ConfigObj(args.inifile)

setParams={}
simParams={}
LRG={'params':{},'wpidx':{},'xi0dx':{},'xi2dx':{},'wp3idx':{},'xi3idx':{},'model':{}}
ELG={'params':{},'wpidx':{},'xi0dx':{},'xi2dx':{},'wp3idx':{},'xi3idx':{},'model':{}}
QSO={'params':{},'wpidx':{},'xi0dx':{},'xi2dx':{},'wp3idx':{},'xi3idx':{},'model':{}}
LXE={'wpidx':{},'xi0idx':{},'xi2idx':{}}
LXQ={'wpidx':{},'xi0idx':{},'xi2idx':{}}
EXQ={'wpidx':{},'xi0idx':{},'xi2idx':{}}
pathIn={}
binParams={'tab':{},'rppi':{},'smu':{},'triXY':{},'tri3D':{}}

tmp=config['setParams']
tmpkey=tmp.keys()
for i in range(len(tmpkey)):
    setParams[tmpkey[i]]=tmp.as_bool(tmpkey[i])
print(setParams)

tmp=config['simParams']
tmpkey=tmp.keys()
for i in range(2):
    simParams[tmpkey[i]]=tmp.as_float(tmpkey[i])
print(simParams)

if(setParams['useLRG']):
    tmp=config['LRG']['params']
    tmpkey=tmp.keys()
    for i in range(len(tmpkey)):
        tmp2=tmp[tmpkey[i]]
        tmpkey2=tmp2.keys()
        LRG['params'][tmpkey[i]]={}
        for j in range(len(tmpkey2)):
            LRG['params'][tmpkey[i]][tmpkey2[j]]=tmp2.as_float(tmpkey2[j])

    tmp=config['LRG']['model']
    tmpkey=tmp.keys()
    for i in range(len(tmpkey)):
        LRG['model'][tmpkey[i]]=tmp.as_int(tmpkey[i])

    if(setParams['useWp']):
        tmp=config['LRG']['wpidx']
        tmpkey=tmp.keys()
        for i in range(len(tmpkey)):
            LRG['wpidx'][tmpkey[i]]=tmp.as_int(tmpkey[i])
    if(setParams['useXil']):
        tmp=config['LRG']['xiidx']
        tmpkey=tmp.keys()
        for i in range(len(tmpkey)):
            LRG['xiidx'][tmpkey[i]]=tmp.as_int(tmpkey[i])
    if(setParams['useWp3']):
        tmp=config['LRG']['wp3idx']
        tmpkey=tmp.keys()
        for i in range(len(tmpkey)):
            LRG['wp3idx'][tmpkey[i]]=tmp.as_int(tmpkey[i])
    if(setParams['useXi3']):
        tmp=config['LRG']['xi3idx']
        tmpkey=tmp.keys()
        for i in range(len(tmpkey)):
            LRG['xi3idx'][tmpkey[i]]=tmp.as_int(tmpkey[i])


    print(LRG)

if(setParams['useELG']):
    tmp=config['ELG']['params']
    tmpkey=tmp.keys()
    for i in range(len(tmpkey)):
        tmp2=tmp[tmpkey[i]]
        tmpkey2=tmp2.keys()
        ELG['params'][tmpkey[i]]={}
        for j in range(len(tmpkey2)):
            ELG['params'][tmpkey[i]][tmpkey2[j]]=tmp2.as_float(tmpkey2[j])
   
    tmp=config['ELG']['model']
    tmpkey=tmp.keys()
    for i in range(len(tmpkey)):
        ELG['model'][tmpkey[i]]=tmp.as_int(tmpkey[i])
    
    if(setParams['useWp']):
        tmp=config['ELG']['wpidx']
        tmpkey=tmp.keys()
        for i in range(len(tmpkey)):
            ELG['wpidx'][tmpkey[i]]=tmp.as_int(tmpkey[i])
    if(setParams['useXil']):
        tmp=config['ELG']['xiidx']
        tmpkey=tmp.keys()
        for i in range(len(tmpkey)):
            ELG['xiidx'][tmpkey[i]]=tmp.as_int(tmpkey[i])
    if(setParams['useWp3']):
        tmp=config['ELG']['wp3idx']
        tmpkey=tmp.keys()
        for i in range(len(tmpkey)):
            ELG['wp3idx'][tmpkey[i]]=tmp.as_int(tmpkey[i])
    if(setParams['useXi3']):
        tmp=config['ELG']['xi3idx']
        tmpkey=tmp.keys()
        for i in range(len(tmpkey)):
            ELG['xi3idx'][tmpkey[i]]=tmp.as_int(tmpkey[i])

   
    print(ELG)

if(setParams['useQSO']):
    tmp=config['QSO']['params']
    tmpkey=tmp.keys()
    for i in range(len(tmpkey)):
        tmp2=tmp[tmpkey[i]]
        tmpkey2=tmp2.keys()
        QSO['params'][tmpkey[i]]={}
        for j in range(len(tmpkey2)):
            QSO['params'][tmpkey[i]][tmpkey2[j]]=tmp2.as_float(tmpkey2[j])

    tmp=config['QSO']['model']
    tmpkey=tmp.keys()
    for i in range(len(tmpkey)):
        QSO['model'][tmpkey[i]]=tmp.as_int(tmpkey[i])
    
    if(setParams['useWp']):
        tmp=config['QSO']['wpidx']
        tmpkey=tmp.keys()
        for i in range(len(tmpkey)):
            QSO['wpidx'][tmpkey[i]]=tmp.as_int(tmpkey[i])
    if(setParams['useXil']):
        tmp=config['QSO']['xiidx']
        tmpkey=tmp.keys()
        for i in range(len(tmpkey)):
            QSO['xiidx'][tmpkey[i]]=tmp.as_int(tmpkey[i])
    if(setParams['useWp3']):
        tmp=config['QSO']['wp3idx']
        tmpkey=tmp.keys()
        for i in range(len(tmpkey)):
            QSO['wp3idx'][tmpkey[i]]=tmp.as_int(tmpkey[i])
    if(setParams['useXi3']):
        tmp=config['QSO']['xi3idx']
        tmpkey=tmp.keys()
        for i in range(len(tmpkey)):
            QSO['xi3idx'][tmpkey[i]]=tmp.as_int(tmpkey[i])

   
    print(QSO)

if(setParams['useLXE']):
    if(setParams['useWp']):
        tmp=config['LXE']['wpidx']
        tmpkey=tmp.keys()
        for i in range(len(tmpkey)):
            LXE['wpidx'][tmpkey[i]]=tmp.as_int(tmpkey[i])
    if(setParams['useXil']):
        tmp=config['LXE']['xiidx']
        tmpkey=tmp.keys()
        for i in range(len(tmpkey)):
            LXE['xiidx'][tmpkey[i]]=tmp.as_int(tmpkey[i])
    print(LXE)
if(setParams['useLXQ']):
    if(setParams['useWp']):
        tmp=config['LXQ']['wpidx']
        tmpkey=tmp.keys()
        for i in range(len(tmpkey)):
            LXQ['wpidx'][tmpkey[i]]=tmp.as_int(tmpkey[i])
    if(setParams['useXil']):
        tmp=config['LXQ']['xiidx']
        tmpkey=tmp.keys()
        for i in range(len(tmpkey)):
            LXQ['xiidx'][tmpkey[i]]=tmp.as_int(tmpkey[i])
    print(LXQ)
if(setParams['useEXQ']):
    if(setParams['useWp']):
        tmp=config['EXQ']['wpidx']
        tmpkey=tmp.keys()
        for i in range(len(tmpkey)):
            EXQ['wpidx'][tmpkey[i]]=tmp.as_int(tmpkey[i])
    if(setParams['useXil']):
        tmp=config['EXQ']['xiidx']
        tmpkey=tmp.keys()
        for i in range(len(tmpkey)):
            EXQ['xiidx'][tmpkey[i]]=tmp.as_int(tmpkey[i])
    print(EXQ)


tmp=config['pathIn']
tmpkey=tmp.keys()
for i in range(len(tmpkey)):
    pathIn[tmpkey[i]]=tmp[tmpkey[i]]
print(pathIn)

tmp=config['binParams']['tab']
tmpkey=tmp.keys()
for i in range(len(tmpkey)):
    binParams['tab'][tmpkey[i]]=tmp.as_int(tmpkey[i])
if(setParams['useWp']):
    tmp=config['binParams']['rppi']
    binParams['rppi']['nrpbins']=tmp.as_int('nrpbins')
    binParams['rppi']['npibins']=tmp.as_int('npibins')
    binParams['rppi']['rpmin']=tmp.as_float('rpmin')
    binParams['rppi']['rpmax']=tmp.as_float('rpmax')
if(setParams['useXil']):
    tmp=config['binParams']['smu']
    tmpkey=tmp.keys()
    for i in range(len(tmpkey)):
        binParams['smu'][tmpkey[i]]=tmp.as_int(tmpkey[i])
if(setParams['useWp3']):
    tmp=config['binParams']['triXY']
    tmpkey=tmp.keys()
    for i in range(len(tmpkey)):
        binParams['triXY'][tmpkey[i]]=tmp.as_int(tmpkey[i])
if(setParams['useXi3']):
    tmp=config['binParams']['tri3D']
    tmpkey=tmp.keys()
    for i in range(len(tmpkey)):
        binParams['tri3D'][tmpkey[i]]=tmp.as_int(tmpkey[i])
 
print(binParams)







