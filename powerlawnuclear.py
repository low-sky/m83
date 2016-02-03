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

idxnuc=np.where(t['RADIUS_PC']<1000)

massnuc=t['MASS_EXTRAP'][idxnuc].data
myfit=powerlaw.Fit(massnuc)

#myfit.plot_ccdf() #look at the data with the fit
print(myfit.alpha) #returns the fitted index

R,p=myfit.distribution_compare('power_law', 'truncated_power_law')
print(R,p)

myfit.truncated_power_law.plot_ccdf(label='Trunc. Power Law')
myfit.power_law.plot_ccdf(label='Power Law')
myfit.plot_ccdf(drawstyle='steps', label='data')
plt.xlabel(r'$Mass\ M_{\odot}$')
plt.ylabel(r'$N\ (>M)$')
plt.legend()
plt.savefig('powerlawnuclear.png')


