from configobj import ConfigObj
config=ConfigObj()
config.filename='example.ini'
config['setParams']={'useLRG':'True',
                     'useELG':'True',
                     'useQSO':'True',
                     'useLXE':'True',
                     'useLXQ':'True',
                     'useEXQ':'True',
                     'useWp':'True',
                     'useXil':'False',
                     'useWp3':'False',
                     'useXi3':'False'}
config['simParams']={'boxsize':'2000.0',
                     'Mpart':'2109081520.453063'}
config['LRG']={
        'params':{
            'logMcut':{'min':'11.0','max':'15.0','init':'13.0','width':'0.01'},
            'logM1':{'min':'11.0','max':'15.0','init':'13.0','width':'0.1'},
            'sigma':{'min':'0.0','max':'2.0','init':'1','width':'0.1'},
            'alpha':{'min':'0.0','max':'2.0','init':'1','width':'0.1'},
            'kappa':{'min':'0.0','max':'5.0','init':'1','width':'0.1'},
            'pmax':{'min':'0.0','max':'1.0','init':'0.3','width':'0.01'}},
        'wpidx':{'min':'4','max':'30'},
        'xi0idx':{'min':'4','max':'30'},
        'xi2idx':{'min':'4','max':'30'},
        'wp3idx':{'min':'4','max':'30'},
        'xi3idx':{'min':'4','max':'30'},
        'model':{'cent':'1','sate':'1'}
       }
config['ELG']={
        'params':{
            'logMcut':{'min':'11.0','max':'15.0','init':'13.0','width':'0.01'},
            'logM1':{'min':'11.0','max':'15.0','init':'13.0','width':'0.1'},
            'sigma':{'min':'0.0','max':'2.0','init':'1','width':'0.1'},
            'alpha':{'min':'0.0','max':'2.0','init':'1','width':'0.1'},
            'kappa':{'min':'0.0','max':'5.0','init':'1','width':'0.1'},
            'pmax':{'min':'0.0','max':'1.0','init':'0.3','width':'0.01'},
            'gamma':{'min':'0.0','max':'10.0','init':'4.0','width':'0.1'},
            'Q':{'min':'0','max':'200','init':'100','width':'1'}},
        'wpidx':{'min':'4','max':'30'},
        'xi0idx':{'min':'4','max':'30'},
        'xi2idx':{'min':'4','max':'30'},
        'wp3idx':{'min':'4','max':'30'},
        'xi3idx':{'min':'4','max':'30'},
        'model':{'cent':'2','sate':'1'}
       }
config['QSO']={
        'params':{
            'logMcut':{'min':'11.0','max':'15.0','init':'13.0','width':'0.01'},
            'logM1':{'min':'11.0','max':'15.0','init':'13.0','width':'0.1'},
            'sigma':{'min':'0.0','max':'2.0','init':'1','width':'0.1'},
            'alpha':{'min':'0.0','max':'2.0','init':'1','width':'0.1'},
            'kappa':{'min':'0.0','max':'5.0','init':'1','width':'0.1'},
            'logMmin':{'min':'11.0','max':'15.0','init':'13.0','width':'0.1'},
            'pmax':{'min':'0.0','max':'1.0','init':'0.3','width':'0.01'}},
        'wpidx':{'min':'13','max':'30'},
        'xi0idx':{'min':'13','max':'30'},
        'xi2idx':{'min':'13','max':'30'},
        'wp3idx':{'min':'13','max':'30'},
        'xi3idx':{'min':'13','max':'30'},
        'model':{'cent':'1','sate':'2'}
        }

config['LXE']={
        'wpidx':{'min':'5','max':'30'}
#        'xi0idx':{'min':'5','max':'30'},
#        'xi2idx':{'min':'5','max':'30'}
        }
config['LXQ']={
        'wpidx':{'min':'5','max':'30'}
#        'xi0idx':{'min':'5','max':'30'},
#        'xi2idx':{'min':'5','max':'30'}
        }
config['EXQ']={
        'wpidx':{'min':'5','max':'30'}
#        'xi0idx':{'min':'5','max':'30'},
#        'xi2idx':{'min':'5','max':'30'}
        }


config['pathIn']={
        'tab':'/mnt/Data2/TabNpcfData_cleaned/rppi/base_c000_ph006/z0.500/halo_mc11_step_part_5p/',
        'rppi':'/mnt/Data2/TabNpcfData_cleaned/rppi/base_c000_ph006/z0.500/halo_mc11_step_part_5p/',
        'smu':'/mnt/Data2/TabNpcfData_cleaned/smu/base_c000_ph006/z0.500/halo_mc11_step_part_5p',
        'triXY':'/mnt/Data2/TabNpcfData_cleaned/triXY/base_c000_ph006/z0.500/halo_mc11_step_part_5p',
        'tri3D':'/mnt/Data2/TabNpcfData_cleaned/tri3D/base_c000_ph006/z0.500/halo_mc11_step_part_5p',
        'wpdata':'wpdata/wp.dat',
        'xildata':'xildata/xil.dat',
        'wp3data':'wp3data/wp3.dat',
        'xi3data':'xil3data/xi3.dat',
        'cov':'data/cov.dat'
        }

config['binParams']={
        'tab':{'nmassbins':20}
        'rppi':{'nrpbins':'25','npibins':'40','rpmin':'-1.5','rpmax':'1.477'},
        'smu':{'nsbins':'25','nmubins':'40'},
        'triXY':{'nsbins':'125'},
        'tri3D':{'nsbins':'125'}
        }

config.write()
