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
import sys

if sys.platform == 'darwin':
    mpl.rcParams['font.family'] = 'serif'
    mpl.rcParams['font.serif'] = 'Times New Roman'

elif sys.platform == 'linux2':
    mpl.rcParams['font.family'] = 'serif'
    mpl.rcParams['font.serif'] = 'FreeSerif'

mpl.rcParams['font.size'] = 12
mygalaxy=Galaxy("M83")

#load the file
mytable=Table.read('m83.co10.K_props_clfind.fits')
#note that matplotlib understands LaTeX math in the code m

rgal=mygalaxy.radius(ra=(mytable['XPOS']), dec=(mytable['YPOS']))
colrgal=Column(name='RADIUS_KPC',data=(rgal/1e3))
mytable.add_column(colrgal)

#idx = np.where(np.isfinite(mytable['VIRMASS_GCORR_DECONV'])) #removes values that are infinite/NaN

#condition = (np.log10(mytable['VIRMASS_GCORR_DECONV']/
#                      mytable['MASS_GCORR'])>-1)
#condition = mytable['PEAK_TO_EDGE']>0
condition = np.ones(len(mytable),dtype=np.bool)
condition = (mytable['MAXVAL']/mytable['NOISE']>5.5)*(rgal>1*u.kpc)*(mytable['RADRMS_GCORR_DECONV']>20)
condition2 = (mytable['MAXVAL']/mytable['NOISE']>5.5)*(rgal<1*u.kpc)*(mytable['RADRMS_GCORR_DECONV']>20)
idx = np.where(condition)
idx2 = np.where(condition2)
mytable2 = mytable[idx2]
mytable = mytable[idx]
#m,b=np.polyfit(np.log(mytable['MASS_GCORR']),np.log(mytable['VIRMASS_GCORR_DECONV']),1)
figure=plt.figure(figsize=(8.0,7.5)) #figure size in inche
ax = plt.subplot(221)
plt.scatter(mytable['MASS_GCORR'],mytable['VIRMASS_GCORR_DECONV'],
            c=mytable['RADIUS_KPC'], marker='s',cmap='Greys',vmin=0,vmax=5)
plt.scatter(mytable2['MASS_GCORR'],mytable2['VIRMASS_GCORR_DECONV'],
            c=mytable2['RADIUS_KPC'], marker='D',cmap='Greys',vmin=0,vmax=5,s=8)
test_mass = np.logspace(5,9,100)
#plt.plot(test_mass,np.exp(m*np.log(test_mass)+b)) #fit line for m83
plt.plot(test_mass,test_mass,alpha=0.5,linewidth=2) #fitted x=y
ax.fill_between([1e5,3.35e5],1e5,5e7, color='gray',alpha=0.5)
virlim = 1040 * 2.57**2 * 20 *4 
ax.fill_between([1e5,5e7],1e5,virlim, color='gray',alpha=0.5)
plt.text(0.05,0.95,'(a)',transform=ax.transAxes,va='top')
plt.xlabel(r'$M_{\mathrm{lum}}\ (M_{\odot})$')
plt.ylabel(r'$M_{\mathrm{vir}}\ (M_{\odot})$')
plt.axis([1e5, 5e7, 1e5, 5e7])
#plt.colorbar(label='$R_{gal}\ (pc)$')
plt.xscale('log')
plt.yscale('log')
plt.tight_layout()
#plt.savefig('MlumMvir_matplotlib.pdf')

#n,c=np.polyfit(np.log(mytable['RADRMS_GCORR_DECONV']),np.log(mytable['VRMS_GCORR_DECONV']),1)
#figure=plt.figure(figsize=(4.5,4)) #figure size in inches
ax2 = plt.subplot(222)
ax2.scatter(mytable['RADRMS_GCORR_DECONV'],mytable['VRMS_GCORR_DECONV'],
            c=mytable['RADIUS_KPC'], marker='s',cmap='Greys',vmin=0,vmax=5)
ax2.scatter(mytable2['RADRMS_GCORR_DECONV'],mytable2['VRMS_GCORR_DECONV'],
            c=mytable2['RADIUS_KPC'], marker='D',cmap='Greys',vmin=0,vmax=5,s=8)


test_radius=np.logspace(0,3,100)
#plt.plot(test_radius, np.exp(n*np.log(test_radius)+c)) #fit line for m83
ax2.plot(test_radius, ((np.pi**0.5/3.4)**0.5)*(test_radius**0.5)) #fit line from Solomon et al.- using effective radius
ax2.text(0.05,0.95,'(b)',transform=ax2.transAxes,va='top')
ax2.set_xlabel(r'$R\ (\mathrm{pc})$')
ax2.set_ylabel(r'${\sigma_v}\  (\mathrm{km}\ \mathrm{s}^{-1})$')
ax2.axis([10,2e2, 1, 3*10e0])
ax2.fill_between([10.,22.8],1e0,30,color='gray',alpha=0.5)
ax2.fill_between([10.,200.],1e0, 2.57,color='gray',alpha=0.5)
plt.colorbar(label='$R_{\mathrm{gal}}\ (\mathrm{kpc})$')
plt.xscale('log')
plt.yscale('log')
plt.tight_layout()
#plt.savefig('Rlinewidth_matplotlib.pdf')

#p,d=np.polyfit(np.log(mytable['RADRMS_GCORR_DECONV']),np.log(mytable['MASS_GCORR']),1)
#figure=plt.figure(figsize=(4.5,4)) #figure size in inches
ax3 = plt.subplot(223)
plt.scatter(mytable['RADRMS_GCORR_DECONV'],mytable['MASS_GCORR'],
            c=mytable['RADIUS_KPC'], marker='s',cmap='Greys',vmin=0,vmax=5)
plt.scatter(mytable2['RADRMS_GCORR_DECONV'],mytable2['MASS_GCORR'],
            c=mytable2['RADIUS_KPC'], marker='D',cmap='Greys',vmin=0,vmax=5,s=8)
test_mr=np.logspace(0,3,100)
#plt.plot(test_mr,np.exp(p*np.log(test_mr)+d)) #fit line for m83
plt.plot(test_mr, (540*test_mr**(2))) #fit line from Solomon et al. - using effective radius
plt.text(0.05,0.95,'(c)',transform=ax3.transAxes,va='top')
plt.ylabel(r'$M_{\mathrm{lum}}\ (M_{\odot})$')
plt.xlabel(r'$R\ (\mathrm{pc})$')
plt.axis([10, 2e2, 10e4, 5e7])
ax3.fill_between([10,22.8],1e5,5e7,color='gray',alpha=0.5)
ax3.fill_between([10,200],1e5,3.35e5,color='gray',alpha=0.5)
#plt.colorbar(label='$R_{gal}\ (pc)$')
plt.xscale('log')
plt.yscale('log')
plt.tight_layout()

ax4 = plt.subplot(224)
plt.scatter((mytable['MASS_GCORR'])/(mytable['RADRMS_GCORR_DECONV'])**2/np.pi,
            (mytable['VRMS_GCORR_DECONV']**2)/(mytable['RADRMS_GCORR_DECONV']),
            marker='s',cmap='Greys',c=mytable['RADIUS_KPC'],vmin=0,vmax=5)
plt.scatter((mytable2['MASS_GCORR'])/(mytable2['RADRMS_GCORR_DECONV'])**2/np.pi,
            (mytable2['VRMS_GCORR_DECONV']**2)/(mytable2['RADRMS_GCORR_DECONV']),
            marker='D',cmap='Greys',c=mytable2['RADIUS_KPC'],vmin=0,vmax=5,s=8)
test_mr=np.logspace(1,4,10)
#plt.plot(test_mr,np.exp(p*np.log(test_mr)+d)) #fit line for m83
plt.plot(test_mr, test_mr/300) #fit line from Solomon et al. - using effective radius
plt.text(0.05,0.95,'(d)',transform=ax4.transAxes,va='top')
plt.ylabel(r'$\sigma_v^2/R\ (\mathrm{km^2\ s^2\ pc}^{-1})$')
plt.xlabel(r'$\Sigma\ (M_{\odot}\ \mathrm{pc}^{-2})$')
plt.colorbar(label='$R_{\mathrm{gal}}\ (\mathrm{kpc})$')
plt.xscale('log')
plt.yscale('log')
plt.ylim(1e-1,1e1)
plt.xlim(1e1,3e3)
plt.tight_layout()


plt.savefig('LarsonLaw.pdf')
