import numpy as np
from astropy.table import Table
from astropy.table import Column
# import matplotlib.pyplot as plt
import powerlaw
from galaxies import Galaxy
import matplotlib as mpl
from sklearn.cluster import k_means

mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = 'Times New Roman'
mpl.rcParams['font.size'] = 12

# get info about m83
mygalaxy = Galaxy("M83")
print(mygalaxy)

# load fits file
tin = Table.read('m83.co10.K_props_clfind.fits')

minmass_global = 4e5
minmass = 4e5
# find cloud's galactocentric distance
tin = tin[tin['MASS_GCORR'] > minmass_global]
rgal = mygalaxy.radius(ra=(tin['XPOS']), dec=(tin['YPOS']))
xgal, ygal = mygalaxy.radius(ra=(tin['XPOS']), dec=(tin['YPOS']),
                             returnXY=True)
_, labels, _ = k_means(np.c_[xgal.value, ygal.value], 20)

colrgal = Column(name='RGAL_PC', data=(rgal))
tin.add_column(colrgal)

t = Table(names=['Xavg', 'Yavg', 'Rmean', 'R_tpl', 'p_tpl', 'index',
                 'index_tpl', 'Mtrunc_tpl',
                 'M1', 'M5', 'Mean5', 'Mmean'])

idx = np.argsort(tin['RGAL_PC'])
tfit = tin[idx]
npts = len(tfit)
stride = 25
mass = tfit['MASS_GCORR']
rgvals = tfit['RGAL_PC']
step = 25

for thislabel in np.unique(labels):
    rads = rgvals[labels == thislabel]
    fitmass = np.sort(mass[labels == thislabel])[::-1]
    minmass = np.max([minmass_global, fitmass.min()])
    fit = powerlaw.Fit(fitmass, xmin=minmass, discrete=True)
    R1, p1 = fit.distribution_compare('power_law', 'truncated_power_law')
    t.add_row()
    t[-1]['Xavg'] = np.mean(tfit[labels == thislabel]['XPOS'])
    t[-1]['Yavg'] = np.mean(tfit[labels == thislabel]['YPOS'])
    t[-1]['Rmean'] = np.mean(rads)
    t[-1]['R_tpl'] = R1
    t[-1]['p_tpl'] = p1
    t[-1]['index'] = fit.alpha
    t[-1]['index_tpl'] = fit.truncated_power_law.parameter1
    t[-1]['Mtrunc_tpl'] = 1 / fit.truncated_power_law.parameter2
    t[-1]['M1'] = fitmass[0]
    t[-1]['M5'] = fitmass[4]
    t[-1]['Mean5'] = np.exp(np.mean(np.log(fitmass[0:5])))
    t[-1]['Mmean'] = np.exp(np.mean(np.log(fitmass)))
t.write('mspecmap_plfit.fits', overwrite=True)

