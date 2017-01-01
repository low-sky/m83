# import libraries
import numpy as np
from astropy.table import Table
from astropy.table import Column
from galaxies import Galaxy
import matplotlib as mpl
import plfit_mcstat as plf
import scipy.stats as ss

mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = 'Times New Roman'
mpl.rcParams['font.size'] = 12


# get info about m83
mygalaxy = Galaxy("M83")
print(mygalaxy)

# load fits file
tin = Table.read('m83.co10.K_props_clfind.fits')

minmass_global = 3.35e5
minmass = 3.35e5
# find cloud's galactocentric distance
rgal = mygalaxy.radius(ra=(tin['XPOS']), dec=(tin['YPOS']))
colrgal = Column(name='RGAL_PC', data=(rgal))
tin.add_column(colrgal)

t = Table(names=['Rmin', 'Rmax', 'Rmean', 'R_tpl', 'p_tpl', 'index',
                 'index_tpl', 'Mtrunc_tpl',
                 'M1', 'M5', 'Mean5', 'Mmean',
                 'index_trunc', 'Mtrun_trunc',
                 'index15_trunc', 'index85_trunc',
                 'Mtrun15_trunc', 'Mtrun85_trunc',
                 'pkprob_trunc'])

tfit = tin[tin['MASS_GCORR'] > minmass_global]
idx = np.argsort(tfit['RGAL_PC'])
tfit = tfit[idx]
npts = len(tfit)
stride = 50
mass = tfit['MASS_GCORR']
rgvals = tfit['RGAL_PC']
step = 25

for i in np.arange(stride, npts - stride, step):
    print i
    fitmass = np.sort(mass[(i - stride):(i + stride)])[::-1]
    rads = rgvals[(i - stride):(i + stride)]
    type = 'trunc'
    optresult = plf.plfit_ksstat(fitmass / minmass, type=type, method='ks')
    sampler = plf.plfit_emcee(fitmass / minmass, type=type, method='ks',
                              steps=200, thin=5)
    pvec = optresult.x
    t.add_row()
    t[-1]['index_' + type] = pvec[0]
    t[-1]['index15_' + type] = ss.scoreatpercentile(
        sampler.flatchain[:, 0], 15)
    t[-1]['index85_' + type] = ss.scoreatpercentile(
        sampler.flatchain[:, 0], 85)
    t[-1]['Mtrun_' + type] = 1e1**pvec[1] * minmass
    t[-1]['Mtrun15_' + type] = 1e1**(ss.scoreatpercentile(
            sampler.flatchain[:, 1], 15)
                                     ) * minmass
    t[-1]['Mtrun85_' + type] = 1e1**(ss.scoreatpercentile(
            sampler.flatchain[:, 1], 85)
                             ) * minmass
    t[-1]['pkprob_' + type] = plf.logprob_plfit(optresult.x,
                                                mass / minmass,
                                                type=type, method='ad')

    t[-1]['Rmin'] = rads.min()
    t[-1]['Rmax'] = rads.max()
    t[-1]['Rmean'] = np.mean(rads)
    t[-1]['M1'] = fitmass[0]
    t[-1]['M5'] = fitmass[4]
    t[-1]['Mean5'] = np.exp(np.mean(np.log(fitmass[0:5])))
    t[-1]['Mmean'] = np.exp(np.mean(np.log(fitmass)))
t.write('rolling_plfit.fits', overwrite=True)
