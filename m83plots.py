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
import matplotlib as mpl
mpl.rcParams['font.family']='serif'
mpl.rcParams['font.serif'] = 'Times New Roman'
mygalaxy=Galaxy("M83")

#load the file
mytable=Table.read('m83.co10.K_props_clfind.fits')
#note that matplotlib understands LaTeX math in the code m

rgal=mygalaxy.radius(ra=(mytable['XPOS']), dec=(mytable['YPOS']))
colrgal=Column(name='RADIUS_KPC',data=(rgal/1e3))
mytable.add_column(colrgal)

#idx = np.where(np.isfinite(mytable['VIRMASS_EXTRAP_DECONV'])) #removes values that are infinite/NaN

#condition = (np.log10(mytable['VIRMASS_EXTRAP_DECONV']/
#                      mytable['MASS_EXTRAP'])>-1)
#condition = mytable['PEAK_TO_EDGE']>0
condition = np.ones(len(mytable),dtype=np.bool)
condition = (mytable['MAXVAL']/mytable['NOISE']>5)*(rgal>1*u.kpc)*(mytable['RADRMS_EXTRAP_DECONV']>20)
condition2 = (mytable['MAXVAL']/mytable['NOISE']>5)*(rgal<1*u.kpc)*(mytable['RADRMS_EXTRAP_DECONV']>20)
idx = np.where(condition)
idx2 = np.where(condition2)
mytable2 = mytable[idx2]
mytable = mytable[idx]
#m,b=np.polyfit(np.log(mytable['MASS_EXTRAP']),np.log(mytable['VIRMASS_EXTRAP_DECONV']),1)
figure=plt.figure(figsize=(8.0,7.5)) #figure size in inche
ax = plt.subplot(221)
plt.scatter(mytable['MASS_EXTRAP'],mytable['VIRMASS_EXTRAP_DECONV'], 
            c=mytable['RADIUS_KPC'], marker='s',cmap='Greys',vmin=0,vmax=5)
plt.scatter(mytable2['MASS_EXTRAP'],mytable2['VIRMASS_EXTRAP_DECONV'], 
            c=mytable2['RADIUS_KPC'], marker='D',cmap='Greys',vmin=0,vmax=5,s=8)
test_mass = np.logspace(5,9,100)
#plt.plot(test_mass,np.exp(m*np.log(test_mass)+b)) #fit line for m83
plt.plot(test_mass,test_mass,alpha=0.5,linewidth=2) #fitted x=y
plt.text(0.05,0.95,'(a)',transform=ax.transAxes,va='top')
plt.xlabel(r'$M_{\mathrm{lum}}\ (M_{\odot})$')
plt.ylabel(r'$M_{\mathrm{vir}}\ (M_{\odot})$')
plt.axis([10e4, 5e7, 10e3, 5e7])
#plt.colorbar(label='$R_{gal}\ (pc)$')
plt.xscale('log')
plt.yscale('log')
plt.tight_layout()
#plt.savefig('MlumMvir_matplotlib.pdf')

#n,c=np.polyfit(np.log(mytable['RADRMS_EXTRAP_DECONV']),np.log(mytable['VRMS_EXTRAP_DECONV']),1)
#figure=plt.figure(figsize=(4.5,4)) #figure size in inches
ax = plt.subplot(222)
plt.scatter(mytable['RADRMS_EXTRAP_DECONV'],mytable['VRMS_EXTRAP_DECONV'], 
            c=mytable['RADIUS_KPC'], marker='s',cmap='Greys',vmin=0,vmax=5)
plt.scatter(mytable2['RADRMS_EXTRAP_DECONV'],mytable2['VRMS_EXTRAP_DECONV'], 
            c=mytable2['RADIUS_KPC'], marker='D',cmap='Greys',vmin=0,vmax=5,s=8)

test_radius=np.logspace(0,3,100)
#plt.plot(test_radius, np.exp(n*np.log(test_radius)+c)) #fit line for m83
plt.plot(test_radius, ((np.pi**0.5/3.4)**0.5)*(test_radius**0.5)) #fit line from Solomon et al.- using effective radius
plt.text(0.05,0.95,'(b)',transform=ax.transAxes,va='top')
plt.xlabel(r'$R\ (\mathrm{pc})$')
plt.ylabel(r'${\sigma_v}\  (\mathrm{km}\ \mathrm{s}^{-1})$')
plt.axis([10,1*10e1, 1, 3*10e0])
plt.colorbar(label='$R_{\mathrm{gal}}\ (\mathrm{kpc})$')
plt.xscale('log')
plt.yscale('log')
plt.tight_layout()
#plt.savefig('Rlinewidth_matplotlib.pdf')

#p,d=np.polyfit(np.log(mytable['RADRMS_EXTRAP_DECONV']),np.log(mytable['MASS_EXTRAP']),1)
#figure=plt.figure(figsize=(4.5,4)) #figure size in inches
ax = plt.subplot(223)
plt.scatter(mytable['RADRMS_EXTRAP_DECONV'],mytable['MASS_EXTRAP'],
            c=mytable['RADIUS_KPC'], marker='s',cmap='Greys',vmin=0,vmax=5)
plt.scatter(mytable2['RADRMS_EXTRAP_DECONV'],mytable2['MASS_EXTRAP'],
            c=mytable2['RADIUS_KPC'], marker='D',cmap='Greys',vmin=0,vmax=5,s=8)
test_mr=np.logspace(0,3,100)
#plt.plot(test_mr,np.exp(p*np.log(test_mr)+d)) #fit line for m83
plt.plot(test_mr, (540*test_mr**(2))) #fit line from Solomon et al. - using effective radius
plt.text(0.05,0.95,'(c)',transform=ax.transAxes,va='top')
plt.ylabel(r'$M_{\mathrm{lum}}\ (M_{\odot})$')
plt.xlabel(r'$R\ (\mathrm{pc})$')
plt.axis([10, 1*10e1, 10e4, 5e7])
#plt.colorbar(label='$R_{gal}\ (pc)$')
plt.xscale('log')
plt.yscale('log')
plt.tight_layout()

ax = plt.subplot(224)
plt.scatter((mytable['MASS_EXTRAP'])/(mytable['RADRMS_EXTRAP_DECONV'])**2/np.pi,
            (mytable['VRMS_EXTRAP_DECONV']**2)/(mytable['RADRMS_EXTRAP_DECONV']), 
            marker='s',cmap='Greys',c=mytable['RADIUS_KPC'],vmin=0,vmax=5)
plt.scatter((mytable2['MASS_EXTRAP'])/(mytable2['RADRMS_EXTRAP_DECONV'])**2/np.pi,
            (mytable2['VRMS_EXTRAP_DECONV']**2)/(mytable2['RADRMS_EXTRAP_DECONV']),
            marker='D',cmap='Greys',c=mytable2['RADIUS_KPC'],vmin=0,vmax=5,s=8)
test_mr=np.logspace(1,4,10)
#plt.plot(test_mr,np.exp(p*np.log(test_mr)+d)) #fit line for m83
plt.plot(test_mr, test_mr/300) #fit line from Solomon et al. - using effective radius
plt.text(0.05,0.95,'(d)',transform=ax.transAxes,va='top')
plt.ylabel(r'$\sigma_v^2/R\ (\mathrm{km^2\ s^2\ pc}^{-1})$')
plt.xlabel(r'$\Sigma\ (M_{\odot}\ \mathrm{pc}^{-2})$')
plt.colorbar(label='$R_{\mathrm{gal}}\ (\mathrm{kpc})$')
plt.xscale('log')
plt.yscale('log')
plt.ylim(1e-1,1e1)
plt.xlim(3e1,1e4)
plt.tight_layout()


plt.savefig('LarsonLaw.pdf')
