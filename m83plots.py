# -*- coding: utf-8 -*-

#a plot of virial mass vs luminous mass
#a plot of radius vs line width
#a plot of luminous mass vs radius

import astropy
import astropy.table
import numpy as np
from astropy.table import Table, Column
import matplotlib.pyplot as plt
from galaxies import Galaxy
import astropy.units as u

mygalaxy=Galaxy("M83")

#load the file
mytable=Table.read('/home/pafreema/Documents/m83.co10.K_props_cprops.fits')
#note that matplotlib understands LaTeX math in the code m

rgal=mygalaxy.radius(ra=(mytable['XPOS']), dec=(mytable['YPOS']))
colrgal=Column(name='RADIUS_PC',data=(rgal))
mytable.add_column(colrgal)

idx = np.where(np.isfinite(mytable['VIRMASS_EXTRAP_DECONV'])) #removes values that are infinite/NaN

m,b=np.polyfit(np.log(mytable['MASS_EXTRAP'][idx]),np.log(mytable['VIRMASS_EXTRAP_DECONV'][idx]),1)
figure=plt.figure(figsize=(4.5,4)) #figure size in inches
plt.scatter(mytable['MASS_EXTRAP'],mytable['VIRMASS_EXTRAP_DECONV'], c=mytable['RADIUS_PC'], marker='s')
test_mass = np.logspace(5,9,100)
plt.plot(test_mass,np.exp(m*np.log(test_mass)+b)) #fit line for m83
plt.plot(test_mass,test_mass,alpha=0.5,linewidth=2) #fitted x=y
plt.xlabel(r'$M_{\mathrm{lum}}\ (M_{\odot})$')
plt.ylabel(r'$M_{\mathrm{vir}}\ (M_{\odot})$')
plt.axis([10e4, 10e7, 10e3, 10e7])
plt.colorbar(label='$R_{gal}\ (pc)$')
plt.xscale('log')
plt.yscale('log')
plt.tight_layout()
plt.savefig('MlumMvir_matplotlib.pdf')

n,c=np.polyfit(np.log(mytable['RADRMS_EXTRAP_DECONV'][idx]),np.log(mytable['VRMS_EXTRAP_DECONV'][idx]),1)
figure=plt.figure(figsize=(4.5,4)) #figure size in inches
plt.scatter(mytable['RADRMS_EXTRAP_DECONV'],mytable['VRMS_EXTRAP_DECONV'], c=mytable['RADIUS_PC'], marker='s')
test_radius=np.logspace(0,3,100)
plt.plot(test_radius, np.exp(n*np.log(test_radius)+c)) #fit line for m83
plt.plot(test_radius, ((np.pi**0.5/3.4)**0.5)*(test_radius**0.5)) #fit line from Solomon et al.- using effective radius
plt.xlabel(r'$R\ (pc)$')
plt.ylabel(r'${\sigma}\  (km\ s^{-1})$')
plt.axis([6,2*10e1, 3*10e-2, 5*10e0])
plt.colorbar(label='$R_{gal}\ (pc)$')
plt.xscale('log')
plt.yscale('log')
plt.tight_layout()
plt.savefig('Rlinewidth_matplotlib.pdf')

p,d=np.polyfit(np.log(mytable['RADRMS_EXTRAP_DECONV'][idx]),np.log(mytable['MASS_EXTRAP'][idx]),1)
figure=plt.figure(figsize=(4.5,4)) #figure size in inches
plt.scatter(mytable['RADRMS_EXTRAP_DECONV'],mytable['MASS_EXTRAP'],c=mytable['RADIUS_PC'], marker='s')
test_mr=np.logspace(0,3,100)
plt.plot(test_mr,np.exp(p*np.log(test_mr)+d)) #fit line for m83
plt.plot(test_mr, (540*test_mr**(2))) #fit line from Solomon et al. - using effective radius
plt.ylabel(r'$M_{\mathrm{lum}}\ (M_{\odot})$')
plt.xlabel(r'$R\ (pc)$')
plt.axis([6, 2*10e1, 10e4, 10e7])
plt.colorbar(label='$R_{gal}\ (pc)$')
plt.xscale('log')
plt.yscale('log')
plt.tight_layout()
plt.savefig('MassRadius_matplotlib.pdf')
