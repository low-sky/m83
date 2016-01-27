# -*- coding: utf-8 -*-

#import libraries
import astropy
import astropy.table
import numpy as np
from astropy.table import Table
import matplotlib.pyplot as plt
from galaxies import Galaxy
import astropy.units as u

mygalaxy=Galaxy("M83")
print(mygalaxy)

mytable=Table.read('/home/pafreema/Documents/m83.co10_props_cprops.fits')

idx = np.where(np.isfinite(mytable['VIRMASS_EXTRAP_DECONV'])) #removes values that are infinite/NaN

rgal=mygalaxy.radius(ra=(mytable['XPOS'][idx]), dec=(mytable['YPOS'][idx]))
print rgal.to(u.pc) #get radii in parsecs

sigma=(mytable['VRMS_EXTRAP_DECONV'][idx])/np.sqrt((mytable['RADRMS_EXTRAP_DECONV'][idx])**2) #sigma naught
dens=(mytable['MASS_EXTRAP'][idx])/(np.pi*((mytable['RADRMS_EXTRAP_DECONV'][idx])**2)) #mass density

m,b=np.polyfit(dens,sigma,1)
x=np.linspace(0,20000,100)
figure=plt.figure(figsize=(4.5,4))
plt.plot (dens, sigma, marker='s', linestyle='None')
plt.plot(x, m*x+b)
plt.xlabel(r'${\Sigma}(M)\ (M_{\odot}/pc^2)$')
plt.ylabel(r'${\sigma}_0\ $')
plt.tight_layout()
plt.savefig('massdensity_matplotlib.png')

figure=plt.figure(figsize=(4.5,4))
plt.plot(rgal, sigma, marker='s', linestyle='None')
plt.xlabel(r'$R_{gal}\ (pc)$')
plt.ylabel(r'${\sigma}_0$')
plt.tight_layout()
plt.savefig('rgal_matplotlib.png')