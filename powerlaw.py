# -*- coding: utf-8 -*-

import astropy
import astropy.table
import numpy as np
import matplotlib.pyplot as plt
from galaxies import Galaxy
import powerlaw
from astropy.table import Table

t=Table.read('/home/pafreema/Documents/m83.co10_props_cprops.fits') #loads into astropy table
mass=t['MASS_EXTRAP'].data #pulls out the mass variable
myfit=powerlaw.Fit(mass)

#myfit.plot_ccdf() #look at the data with the fit
print(myfit.alpha) #returns the fitted index

R,p=myfit.distribution_compare('power_law', 'truncated_power_law')
print(R,p)
#compares two different forms of the distribution
#R<0 means the data is more consistent with the right distribution, p is the significance of the result
myfit.truncated_power_law.plot_ccdf(label='Trunc. Power Law')
myfit.power_law.plot_ccdf(label='Power Law')
myfit.plot_ccdf(drawstyle='steps', label='data')
plt.xlabel(r'$Mass\ M_{\odot}$')
plt.ylabel(r'$N\ (>M)$')
plt.legend()
plt.savefig('powerlaw.png')
