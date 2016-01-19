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
mytable=Table.read('/home/pafreema/Documents/m83.co10_props_cprops.fits')
#note that matplotlib understands LaTeX math in the code below

m,b=np.polyfit(np.log(mytable['MASS_EXTRAP']),np.log(mytable['VIRMASS_EXTRAP_DECONV']),1)
figure=plt.figure(figsize=(4.5,4)) #figure size in inches
plt.loglog(mytable['MASS_EXTRAP'],mytable['VIRMASS_EXTRAP_DECONV'],marker='s',linestyle='None')
plt.plot(np.log(mytable['MASS_EXTRAP']),m*(np.log(mytable['VIRMASS_EXTRAP_DECONV']))+b,"-r")
plt.xlabel(r'$M_{\mathrm{lum}}\ (M_{\odot})$')
plt.ylabel(r'$M_{\mathrm{vir}}\ (M_{\odot})$')
plt.tight_layout()
plt.savefig('MlumMvir_matplotlib.png')

figure=plt.figure(figsize=(4.5,4)) #figure size in inches
plt.loglog(mytable['RADRMS_EXTRAP_DECONV'],mytable['VRMS_EXTRAP_DECONV'],marker='s',linestyle='None')
plt.xlabel(r'$R\ (pc)$')
plt.ylabel(r'${\sigma}\  (km\ s^{-1})$')
plt.tight_layout()
plt.savefig('Rlinewidth_matplotlib.png')

figure=plt.figure(figsize=(4.5,4)) #figure size in inches
plt.loglog(mytable['MASS_EXTRAP'],mytable['RADRMS_EXTRAP_DECONV'],marker='s',linestyle='None')
plt.xlabel(r'$M_{\mathrm{lum}}\ (M_{\odot})$')
plt.ylabel(r'$R\ (pc)$')
plt.tight_layout()
plt.savefig('MlumR_matplotlib.png')