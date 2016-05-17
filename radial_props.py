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

mytable=Table.read('m83.co10.K_props_clfind.fits')
#note that matplotlib understands LaTeX math in the code m

rgal=mygalaxy.radius(ra=(mytable['XPOS']), dec=(mytable['YPOS']))
colrgal=Column(name='RADIUS_KPC',data=(rgal/1e3))
mytable.add_column(colrgal)

condition = np.ones(len(mytable),dtype=np.bool)
condition = (mytable['MAXVAL']/mytable['NOISE']>5) & (mytable['RADRMS_EXTRAP_DECONV']>20)
idx = np.where(condition)
mytable = mytable[idx]

s2n = np.log10(mytable['MAXVAL']/mytable['NOISE'])

plt.clf()
#plt.figure(figsize=(8.0,4.25))
fig, axes = plt.subplots(nrows=1, ncols=2,figsize=(8,3.5))
fig.subplots_adjust(wspace=0.3,bottom=0.2,right=0.85)
#plt.subplot(211)
ax0 = axes[0]
im1 = ax0.scatter(mytable['RADIUS_KPC'],mytable['VRMS_EXTRAP_DECONV']/mytable['RADRMS_EXTRAP_DECONV']**0.5,marker='s',cmap='Greys',color='grey',c=np.log10(mytable['MASS_EXTRAP']))
ax0.set_yscale('log')
ax0.set_ylabel(r'$\sigma_v/\sqrt{R}\ (\mathrm{km\ s^{-1}\ pc^{-1/2}})$')
ax0.set_xlabel(r'$R_{\mathrm{gal}}\ (\mathrm{kpc})$')
ax0.set_xlim([-0.5,6])

ax1 = axes[1]
im2 = ax1.scatter(mytable['RADIUS_KPC'],mytable['MASS_EXTRAP']/mytable['RADRMS_EXTRAP_DECONV']**2/np.pi,marker='s',cmap='Greys',color='grey',c=np.log10(mytable['MASS_EXTRAP']))
ax1.set_yscale('log')
ax1.set_xlim([-0.5,6])
ax1.set_xlabel(r'$R_{\mathrm{gal}}\ (\mathrm{kpc})$')
ax1.set_ylabel(r'$\Sigma\ (M_{\odot}\ \mathrm{pc}^{-2})$')

#cb =fig.colorbar(im2, ax=axes.ravel().tolist())
cbar_ax =  fig.add_axes([0.86, 0.2, 0.025, 0.7])
cb = fig.colorbar(im2, cax=cbar_ax)
cb.set_label(r'$\log_{10}(M/M_\odot)$')
#fig.tight_layout(pad=5.0)
fig.savefig('RadialPlots.pdf')
