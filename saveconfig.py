from configobj import ConfigObj
config=ConfigObj()
config.filename='example.ini'
config['settings']={'doFitting':'True','HODcatalog':'True'}
config['fit']={
        'params':{
            'logMcut':{'min':'11.0','max':'15.0','init':'13.0'},
            'logM1':{'min':'11.0','max':'15.0','init':'13.0'},
            'sigma':{'min':'0.0','max':'2.0','init':'1'},
            'alpha':{'min':'0.0','max':'2.0','init':'1'},
            'kappa':{'min':'0.0','max':'2.0','init':'1'}},
        'use2pt':{
                'wp2':{
                    'usewp2':'True',
                    'doPaircounts':'True',
                    'inputPaircounts':'False',
                    'logrpmin':'-1.5',
                    'logrpmax':'1.477',
                    'nrpbins':'25',
                    'npibins':'40',
                    'pimax':'40'}},
        'use3pt':{
                'xi3':{
                    'usexi3':'False',
                    'inputTrianglecounts':'False',
                    'nsbins':'125'}},
        'common':{
            'nmassbins':'20',
            'MCMCsteps':'20000'},
        'chainsOut':'chains/',
        'covIn':'cov/',
        'data2In':'data/',
        'data3In':'data/'}

config['readCounts']={
        'massbinmidIn':'massbin/',
        'nummassbinhaloIn':'massbin/',
        'nummassbinpartIn':'massbin/',
        'paircountsIn':'paircountsIn/',
        'trianglecountsIn':'trianglecountsIn/'}

config['paircounts']={
        'common':{
            'nthreads':'8',
            'boxsize':'2000.0',
            'savepc':'True'},
        'paircountsOut':'paircounts/',
        'simulationIn':'/mnt/Data/Abacus/z0.500/halo_info/'}

config['HODcat']={
        'logMcut':'13.0',
        'logM1':'14.0',
        'sigma':'1.0',
        'kappa':'0.5',
        'alpha':'1.0'}

config['cosmology']={'Om':'0.315192','redshift':'0.5','particleMass':'56945201052.2327'}
config.write()
