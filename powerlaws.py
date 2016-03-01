# -*- coding: utf-8 -*-

import astropy
import astropy.table
import numpy as np
import matplotlib.pyplot as plt
from galaxies import Galaxy
import powerlaw
from astropy.table import Table

t=Table.read('/home/pafreema/Documents/m83.co10.K_props_cprops.fits') #loads into astropy table
mass=t['MASS_EXTRAP'].data #pulls out the mass variable
myfit=powerlaw.Fit(mass)
print(myfit.xmin)
idx=np.where(t['MASS_EXTRAP']>3e5)
mass_subset=t['MASS_EXTRAP'][idx].data
myfit_subset=powerlaw.Fit(mass_subset, xmin=3e5)

print(myfit.xmin)

#myfit.plot_ccdf() #look at the data with the fit
print(myfit.alpha) #returns the fitted index
print(myfit_subset.alpha)

R,p=myfit_subset.distribution_compare('power_law', 'truncated_power_law')
print(R,p)
print(1/myfit_subset.truncated_power_law.parameter2)
#compares two different forms of the distribution
#R<0 means the data is more consistent with the right distribution, p is the significance of the result
myfit_subset.truncated_power_law.plot_ccdf(label='Trunc. Power Law')
myfit_subset.power_law.plot_ccdf(label='Power Law')
myfit_subset.plot_ccdf(drawstyle='steps', label='data')
plt.xlabel(r'$Mass\ M_{\odot}$')
plt.ylabel(r'$N\ (>M)$')
plt.legend(loc=3)
plt.savefig('powerlaw.png')