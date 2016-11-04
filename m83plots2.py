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

mytable=Table.read('/home/pafreema/Documents/m83.co10.K_props_clfind.fits')

idx = np.where(np.isfinite(mytable['VIRMASS_GCORR_DECONV'])) #removes values that are infinite/NaN

rgal=mygalaxy.radius(ra=(mytable['XPOS'][idx]), dec=(mytable['YPOS'][idx]))
#print rgal #get radii in parsecs

sigma=(mytable['VRMS_GCORR_DECONV'][idx])/np.sqrt((mytable['RADRMS_GCORR_DECONV'][idx])**2) #sigma naught
dens=(mytable['MASS_GCORR'][idx])/(np.pi*((mytable['RADRMS_GCORR_DECONV'][idx])**2)) #mass density

m,b=np.polyfit(np.log(dens),np.log(sigma),1)
x=np.linspace(10,10e4,100)
figure=plt.figure(figsize=(4.5,4))
plt.scatter(dens, sigma,c=rgal, marker='s')
plt.plot(x, np.exp(m*np.log(x)+b))
plt.xlabel(r'${\Sigma}M\ (M_{\odot}/pc^2)$')
plt.ylabel(r'${\sigma}_0\ $')
plt.axis([3*10, 10e3,10e-3, 1 ])
plt.xscale('log')
plt.yscale('log')
plt.colorbar(label='$R_{gal}\ (pc)$')
plt.tight_layout()
plt.savefig('massdensity_matplotlib.pdf')

figure=plt.figure(figsize=(4.5,4))
plt.plot(rgal, sigma, marker='s', linestyle='None')
plt.xlabel(r'$R_{gal}\ (pc)$')
plt.ylabel(r'${\sigma}_0$')
plt.tight_layout()
plt.savefig('rgal_matplotlib.pdf')
