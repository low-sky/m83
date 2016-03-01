# -*- coding: utf-8 -*-

#a plot of virial mass vs luminous mass
#a plot of radius vs line width
#a plot of luminous mass vs radius

import astropy
import astropy.table
import numpy as np
#import some libraries
from astropy.table import Table
import matplotlib.pyplot as plt

#load the file
mytable=Table.read('/home/pafreema/Documents/m83.co10.K_props_cprops.fits')
#note that matplotlib understands LaTeX math in the code m

idx = np.where(np.isfinite(mytable['VIRMASS_EXTRAP_DECONV'])) #removes values that are infinite/NaN
m,b=np.polyfit(np.log(mytable['MASS_EXTRAP'][idx]),np.log(mytable['VIRMASS_EXTRAP_DECONV'][idx]),1)
figure=plt.figure(figsize=(4.5,4)) #figure size in inches
plt.loglog(mytable['MASS_EXTRAP'],mytable['VIRMASS_EXTRAP_DECONV'],marker='s',linestyle='None')
test_mass = np.logspace(5,8,100)
plt.plot(test_mass,np.exp(m*np.log(test_mass)+b)) #fit line for m83
plt.plot(test_mass,test_mass,alpha=0.5,linewidth=2) #fitted x=y
plt.xlabel(r'$M_{\mathrm{lum}}\ (M_{\odot})$')
plt.ylabel(r'$M_{\mathrm{vir}}\ (M_{\odot})$')
plt.tight_layout()
plt.savefig('MlumMvir_matplotlib.png')

n,c=np.polyfit(np.log(mytable['RADRMS_EXTRAP_DECONV'][idx]),np.log(mytable['VRMS_EXTRAP_DECONV'][idx]),1)
figure=plt.figure(figsize=(4.5,4)) #figure size in inches
plt.loglog(mytable['RADRMS_EXTRAP_DECONV'],mytable['VRMS_EXTRAP_DECONV'],marker='s',linestyle='None')
test_radius=np.logspace(0,3,100)
plt.plot(test_radius, np.exp(n*np.log(test_radius)+c)) #fit line for m83
plt.plot(test_radius, ((np.pi**0.5/3.4)**0.5)*(test_radius**0.5)) #fit line from Solomon et al.- using effective radius
plt.xlabel(r'$R\ (pc)$')
plt.ylabel(r'${\sigma}\  (km\ s^{-1})$')
plt.tight_layout()
plt.savefig('Rlinewidth_matplotlib.png')

p,d=np.polyfit(np.log(mytable['RADRMS_EXTRAP_DECONV'][idx]),np.log(mytable['MASS_EXTRAP'][idx]),1)
figure=plt.figure(figsize=(4.5,4)) #figure size in inches
plt.loglog(mytable['RADRMS_EXTRAP_DECONV'],mytable['MASS_EXTRAP'],marker='s',linestyle='None')
test_mr=np.logspace(0,3,100)
plt.plot(test_mr,np.exp(p*np.log(test_mr)+d)) #fit line for m83
plt.plot(test_mr, (540*test_mr**(2))) #fit line from Solomon et al. - using effective radius
plt.ylabel(r'$M_{\mathrm{lum}}\ (M_{\odot})$')
plt.xlabel(r'$R\ (pc)$')
plt.tight_layout()
plt.savefig('MassRadius_matplotlib.png')
