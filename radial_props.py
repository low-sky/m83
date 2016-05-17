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
plt.figure(figsize=(4.25,5))
plt.subplot(211)
plt.scatter(mytable['RADIUS_KPC'],mytable['VRMS_EXTRAP_DECONV']/mytable['RADRMS_EXTRAP_DECONV']**0.5,marker='s',cmap='Greys',color='grey',c=s2n)
plt.yscale('log')
plt.ylabel(r'$\sigma_v/\sqrt{R}\ (\mathrm{km\ s^{-1}\ pc^{-1/2}})$')
plt.xlim([-0.5,6])
plt.subplot(212)
plt.scatter(mytable['RADIUS_KPC'],mytable['MASS_EXTRAP']/mytable['RADRMS_EXTRAP_DECONV']**2/np.pi,marker='s',cmap='Greys',color='grey',c=s2n)
plt.yscale('log')
plt.xlim([-0.5,6])
plt.xlabel(r'$R_{\mathrm{gal}}\ (\mathrm{kpc})$')
plt.ylabel(r'$\Sigma\ (M_{\odot}\ \mathrm{pc}^{-2})$')
plt.tight_layout()
plt.savefig('RadialPlots.pdf')
