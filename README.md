# TabNpcf-HOD
 
Tabulated N Point Correlation Functions HOD fitting pipeline based on Abacus Simulation

`TabNpcf-HOD` is a fast HOD fitting pipeline using particle based tabulation method. 

Dependencies:

* [Corrfunc](https://corrfunc.readthedocs.io/en/master/index.html)
* [Abacus Utils (compaso_halo_catalog module)](https://abacusutils.readthedocs.io/en/latest/compaso.html#module-abacusnbody.data.compaso_halo_catalog)
* [numba](http://numba.pydata.org/)
* [configobj](https://pypi.org/project/configobj/)
* numpy,scipy,glob


## Quick Start

`TabNPCF-HOD` could be used as a black box for standard HOD fitting. For now only wp part are functional.

Users can take advantage of our pre-computed tabulation paircounts stored on NERSC and available for DESI members. 
We provided precomputed rp-pi paircounts for projected CF based on AbacusSummit_base_c000_ph006 cleaned catalog.
We have them ready for all primary redshifts. The data product are located in the following directory:
> /global/cscratch1/sd/hanyuz/TabNpcfData_cleaned

Paircounts are based on [Corrfunc](https://corrfunc.readthedocs.io/en/master/index.html) library. 
The binning for these paircounts are as follow:

```python
rpbins=numpy.logspace(-1.5,1.477,26)
pibins=numpy.linspace(0,40,41)
```
we have 25 rp bins from around 0.0316 (10^-1.5) to 30 (10^1.477) Mpc/h and maximum separation of 40 Mpc/h along the Z-dimension.

The projected correlation function you are fitting to must use same binning as the pre-computed paircounts. 
However, you can select a range of bins to use for different tracers, so these paircounts are flexible enough for most cases.


Users only need to
* prepare a wp measurements file you are fitting to. (match the binning above, list as order LRG ELG QSO LXE LXQ EXQ, we read second column by default)
* prepare a full cov-matrix file for likelihood calculation.(match the measurements order, if you have 50 points to fit, then the cov-mat should be 50x50)
* modify .ini file for MCMC accordingly (set measurements and cov-mat input path, set bin index,select tracers to use, select HOD model, set HOD params flat prior, init value and estimated error)

Please check INI file for MCMC section and HOD model section for more detail about setting.

```bash
$ python driver.py -ini example_mcmc.ini
```

This is a test 'mock to mock' fitting using baseline HOD model. The MCMC chains and a contour will be saved at:

> chains/example_chains.h5\
> plots/example_plot.pdf

You could then use `getdist` to further analysis the chains.


s-mu paircounts for 2PCF multipole and triangle counts for higher order statistic will be updated later. 

## Full Procedure

For people want to use different simulation box or binning setting, you might check procedure as follow.

### Down sampling

First we need to apply a down sampling to the halo and particle catalog. AbacusSummit have halos and a 3% particles list output.\
However, for tabulation method, we don't need all low-mass halos and all 3% particles for a precise estimation of clustering.

#### INI file for downsampling

```bash
abacusbox=AbacusSummit_base_c000_ph006      #Abacus box name
redshift=z0.500                             #Redshift,support primary redshift only
outpath=downsample/                         #Output directory for down sampled catalog
partfrac=0.05                               #Fraction of particles, partfrac=0.05 will take 0.05*0.03=0.15% of total particles.
randseed=123                                #Randseed for downsample and particle based HOD random list.
masscut=1e11                                #Hard mass cutoff boundary for halo list
halostep=True                               #Considering using step function to futher downsample halo list or not
usersd=True                                 #Considering redshift distortion or not
```

Find an example at config/down.ini. After setting ini file, run as follow to generate downsampled catalog.
```bash
$ export HDF5_USE_FILE_LOCKING=FALSE
$ python prepDown.py -ini config/down.ini
```

### Tabulation and Paircount

Second we tabulate the downsampled catalog and compute Halo-Halo, Halo-Part and Part-Part paircounts cross massbins.

#### INI file for paircounting
```bash
clustering=rppi                             #clustering type, rppi for projected CF, smu for 2PCF multipole (will updated later)
downcatpath=downsample/                     #downsampled catalog directory, the one you set at the first step. 
outpath=rppi/                               #output path for HH HP and PP paircounts
nfiles=9                                    #number of downsampled catalog files, default value is 9.
nrpbins=25                                  #number of rp bins
rpmin=-1.5                                  #lower bound of rp, log10 value
rpmax=1.477                                 #upper bound of rp, log10 value
pimax=40                                    #maximum separation along the Z-dimension
nthreads=32                                 #number of threads for corrfunc library
boxsize=2000.0                              #simulation boxsize
nmassbins=20                                #number of mass bins to tabulate the catalog
particlemass=2109081520.453063              #particle mass of simulation box, default value for base Abacus box are 2109081520.453063
```

Finding an example at config/pair.ini. Then run
```bash
$ python paircount.py -ini config/pair.ini
```

### MCMC
After we get all HH, HP and PP paircounts, we can do MCMC fitting.

#### INI file for MCMC
```bash
[setParams]
useLRG/ELG/QSO/LXE/LXQ/EXQ = True   #Decide with tracer to use.
useWp/Xil/Wp3/Xi3                   #Decide which clustering statistics to use, support WP only for now.

[simParams]
boxsize = 2000.0                    #simulation boxsize
Mpart = 2109081520.453063           #simulation particle mass

[LRG/ELG/QSO/LXE/LXQ/EXQ]
[[params]]
[[[logMcut/logM1/sigma/alpha...(hod parameters, different for different tracer)]]]
min = 11.0                          #params are restricted to stay within the bounds determined by min and max.
max = 15.0                          #params are restricted to stay within the bounds determined by min and max.
init = 13                           #params init value
width = 0.01                        #estimate error for params
[[wpidx]]
min = 4                             #the rp index we want to use for fitting, paircounts might have wide rp range
max = 25                            #we will use only [min:max] points.
[[model]]                           #check HOD model section for detail
cent=1                              #we provide different hod model for different tracer, cent for central HOD
sate=1                              #we provide different hod model for different tracer, cent for satellite HOD

[pathIn]
tab=rppi/                           #numHaloBin.dat etc. file directory. Paircount step will output these files.
rppi=rppi/                          #HH HP PP rppi paircounts directory. Paircount step will output these files.
smu=smu/                            #HH HP PP smu paircounts directory. Paircount step will output these files.
triXY=trixy/                        #HHH HPP PHH PPP triangle counts directory. Provided later by author.
tri3D=tri3D/                        #HHH HPP PHH PPP triangle counts directory. Provided later by author.
wpdata=data/wp.dat                  #wp measurement to fit, we read second column by default. 
                                    please remove all bins you don't want to use, 
                                    we do not check if bins are match to the paircount.
                                    list all wp measurements you want to use as this order:
                                    LRG, ELG, QSO, LXE, LXQ, EXQ.
xildata=data/xil.dat                #xil measurements to fit, save as wp
wp3data=data/wp3.dat                #wp3 measurements to fit, we read fourth column by default
xi3data=data/xi3.dat                #xi3 measurements to fit, we read fourth column by default
cov=data/cov.dat                    #Full cov matrix, ntotalpoints by ntotalpoints.

[binParams]
[[tab]]
nmassbins=20                        #number of mass bins when you tabulate the catalog
[[rppi]]
nrpbins=25                          #rp pi binning setting, should same as paircounts setting.
npibins=40
rpmin=-1.5
rpmax1.477

[mcmc]
steps=20000                         #number of steps for MCMC
outpath                             #path to save MCMC chains
```

Run as follow to do MCMC fitting.
```
$ python driver.py -ini config/example.ini
```

## HOD Model

We include some common HOD model for fitting. Please contact author or modify `hodmodel.py` and `driver.py` to add more.

### Baseline HOD

`Ncent=1.0/2.0*erfc((log10Mcut-log10(Mhalo))/sqrt(2)/log10(e)/sigma)`

`Nsate=((Mhalo-kappa*Mcut)/M1)**alpha`

We have 5 parameters here `log(Mcut), sigma, log(M1), kappa, alpha`

### LRG

Set `[LRG] [[model]] cent=1` in .ini file for MCMC, `Ncent(lrg)=pmax*Ncent`. Additional parameters `pmax` is fixed to 1 for LRG only fitting, is free parameters for multitracer fitting.

Set `[LRG] [[model]] sate=1` in .ini file for MCMC, `Nsate(lrg)=Ncent(lrg)*Nsate`.

### ELG

Set `[ELG] [[model]] cent=1` in .ini file for MCMC, `Ncent(elg)=pmax*Ncent`. `pmax` are free parameter for both single and multitracer case.

Set `[ELG] [[model]] cent=2` in .ini file for MCMC, `Ncent(elg)` follow Shadab et al. [MTHOD](https://arxiv.org/abs/1910.05095) equation 9 as HMQ model. kappa will be fixed to 1 in this case.

Set `[ELG] [[model]] sate=1` in .ini file for MCMC, `Nsate(ELG)=Nsate`.
 
### QSO

Set `[QSO] [[model]] cent=1` in .ini file for MCMC, `Ncent(QSO)=pmax*Ncent`. `pmax` are free parameter for both single and multitracer case.

Set `[QSO] [[model]] sate=1` in .ini file for MCMC, `Nsate(QSO)=Nsate`.

Set `[QSO] [[model]] sate=2` in .ini file for MCMC, `Nsate(QSO)=(Mhalo/M1)**alpha*exp(-Mcut/Mhalo)`.

### Multi-tracer

If you set all sate and cent to 1, it will use Shadab et al. [MTHOD](https://arxiv.org/abs/1910.05095) ELG(erf model) (18 parameters in total), if you change ELG cent to 2, it will use [MTHOD](https://arxiv.org/abs/1910.05095) ELG(HMQ) model (19 parameters in total).

