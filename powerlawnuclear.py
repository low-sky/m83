# -*- coding: utf-8 -*-

#import libraries
import astropy
import astropy.table
import numpy as np
from astropy.table import Table
from astropy.table import Column
import matplotlib.pyplot as plt
import powerlaw
from galaxies import Galaxy
import astropy.units as u

mygalaxy=Galaxy("M83")
print(mygalaxy)

t=Table.read('/home/pafreema/Documents/m83.co10_props_cprops.fits')

rgal=mygalaxy.radius(ra=(t['XPOS']), dec=(t['YPOS']))

colrgal=Column(name='RADIUS_PC',data=(rgal))
t.add_column(colrgal)
#print(t)

#fit for the whole galaxy
mass=t['MASS_EXTRAP'].data #pulls out the mass variable
myfit=powerlaw.Fit(mass)

idxnuc=np.where(t['RADIUS_PC']<1000)
idxmass=np.where((t['MASS_EXTRAP']>2496477.5)&(t['RADIUS_PC']<1000))
massnuc=t['MASS_EXTRAP'][idxnuc].data
massnuc_subset=t['MASS_EXTRAP'][idxmass].data
fit=powerlaw.Fit(massnuc)
fit_subset=powerlaw.Fit(massnuc_subset, xmin=2496477.5)

#myfit.plot_ccdf() #look at the data with the fit
print(fit.alpha) #returns the fitted index
print(fit_subset.alpha)

R,p=fit_subset.distribution_compare('power_law', 'truncated_power_law')
print(R,p)

fit_subset.truncated_power_law.plot_ccdf(label='Trunc. Power Law')
fit_subset.power_law.plot_ccdf(label='Power Law')
fit_subset.plot_ccdf(drawstyle='steps', label='data')
plt.xlabel(r'$Mass\ M_{\odot}$')
plt.ylabel(r'$N\ (>M)$')
plt.legend()
plt.savefig('powerlawnuclear.png')


