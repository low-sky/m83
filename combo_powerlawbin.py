# -*- coding: utf-8 -*-

#import libraries
import astropy
# import astropy.table
import numpy as np
from astropy.table import Table
from astropy.table import Column
import matplotlib.pyplot as plt
# import powerlaw
from galaxies import Galaxy
# import astropy.units as u
import matplotlib as mpl
import astropy.io.fits as fits
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = 'Times New Roman'
mpl.rcParams['font.size'] = 14
#get info about m83
mygalaxy = Galaxy("M83")
print(mygalaxy)

# load fits file
t = Table.read('m83.co10.K_props_clfind.fits')
# import astropy.wcs as wcs

data = fits.getdata('m83.co10.tmax.fits')
hdr = fits.getheader('m83.co10.tmax.fits')
rgal_img = (mygalaxy.radius(header=hdr)).value

# find cloud's galactocentric distance


rgal = mygalaxy.radius(ra=(t['XPOS']), dec=(t['YPOS']))
colrgal = Column(name='RGAL_PC', data=(rgal))
t.add_column(colrgal)

t = t[(t['MASS_EXTRAP'].data > 3e5)]


massvals = t['MASS_EXTRAP'].data
radvals = t['RGAL_PC'].data
radsort = np.sort(radvals)
msort = massvals[np.argsort(radvals)]
mcum = np.cumsum(msort)
massedges = np.linspace(0, 1, 6) * mcum.max()
radedges = np.zeros(massedges.size)
for idx, edge in enumerate(massedges):
    radedges[idx] = np.min(radsort[mcum >= edge])


edges = np.array([0, 450, 2300, 3200, 3900, 4500])
edges = radedges

plt.clf()
plt.figure(figsize=(5.25, 4.5))
apix = (hdr['CDELT2'] * np.pi / 180 * 4.8e6 / 1e3)**2 *\
    np.cos(mygalaxy.inclination)
shapearray = ['o', '+', 's', 'x', 'd']
for ins,outs,shape in zip(edges[0:-1],edges[1:],shapearray):
    subset = (t['MASS_GCORR']>3e5)*(t['RGAL_PC']<=outs)*(t['RGAL_PC']>ins)
    tt = t[subset]
    mass = tt['MASS_GCORR'].data
    area = np.sum((rgal_img>ins)*(rgal_img)<=outs)*apix

    plt.loglog(np.sort(mass),np.linspace(len(mass),1,len(mass))/area,
               label=r'${0:2.1f}<R_g/\rm kpc<{1:2.1f}$'.format(ins/1e3,outs/1e3),marker=shape)
plt.legend(fontsize=9,loc='lower left',numpoints=2)
plt.xlabel(r'$M\ (M_{\odot})$')
plt.ylabel(r'$n(M^\prime>M)\ (\mathrm{kpc}^{-2})$')
plt.ylim(1e-2,2e2)
plt.tight_layout()
plt.savefig('m83_massdist.pdf')

# #add those distances to the fits table
# #print(t)

# #fit for the whole galaxy
# mass=t['MASS_GCORR'].data
# myfit=powerlaw.Fit(mass)
# idx=np.where(t['MASS_GCORR']>3e5)
# mass_subset=t['MASS_GCORR'][idx].data
# myfit_subset=powerlaw.Fit(mass_subset, xmin=3e5)
# R,p=myfit_subset.distribution_compare('power_law', 'truncated_power_law')

# #plot bin 1

# #index the file so you get clouds within a certain radius
# idx_1=np.where(t['RADIUS_PC']<450)
# idxmass_1=np.where((t['MASS_GCORR']>3e5)&(t['RADIUS_PC']<450))

# #pull out the mass variable and fit it
# massdisk1=t['MASS_GCORR'][idx_1].data
# fit_1=powerlaw.Fit(massdisk1)
# massdisk1_subset=t['MASS_GCORR'][idxmass_1].data
# fit_1_subset=powerlaw.Fit(massdisk1_subset, xmin=3e5)


# #compares two different forms of the distribution
# #R<0 means the data is more consistent with the right distribution, p is the significance of the result
# R1,p1=fit_1_subset.distribution_compare('power_law', 'truncated_power_law')

# fig=plt.figure(figsize=(7.5,5))
# ax=fig.add_subplot(111)
# ax.spines['top'].set_color('none')
# ax.spines['bottom'].set_color('none')
# ax.spines['left'].set_color('none')
# ax.spines['right'].set_color('none')
# ax.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')
# ax.set_ylabel('$N\ (>M)$')
# ax.set_xlabel('Mass $(M_{\odot})$')

# ax1=fig.add_subplot(231)
# fit_1_subset.truncated_power_law.plot_ccdf(label='Trunc. Power Law')
# fit_1_subset.power_law.plot_ccdf(label='Power Law')
# fit_1_subset.plot_ccdf(drawstyle='steps', label='Data')
# plt.axis([10e4, 10e7, 10e-4, 10e-1])
# plt.text(2*10e4, 2*10e-4, '(a)')
# plt.setp(ax1.get_xticklabels(), visible=False)
# #plt.legend(loc=3)
# #plt.savefig('powerlawbin1.png')

# #plot bin 2

# idx_2=np.where((t['RADIUS_PC']>450)&(t['RADIUS_PC']<2300))
# idxmass_2=np.where((t['MASS_GCORR']>3e5)&(t['RADIUS_PC']>450)&(t['RADIUS_PC']<2300))

# massdisk2_subset=t['MASS_GCORR'][idxmass_2].data
# fit_2_subset=powerlaw.Fit(massdisk2_subset, xmin=3e5)
# massdisk2=t['MASS_GCORR'][idx_2].data
# fit_2=powerlaw.Fit(massdisk2)

# R2,p2=fit_2_subset.distribution_compare('power_law', 'truncated_power_law')

# ax2=fig.add_subplot(232, sharex=ax1, sharey=ax1)
# fit_2_subset.truncated_power_law.plot_ccdf(label='Trunc. Power Law')
# fit_2_subset.power_law.plot_ccdf(label='Power Law')
# fit_2_subset.plot_ccdf(drawstyle='steps', label='Data')
# plt.axis([10e4, 10e7, 10e-4, 10e-1])
# plt.text(2*10e4, 2*10e-4, '(b)')
# plt.setp(ax2.get_xticklabels(), visible=False)
# plt.setp(ax2.get_yticklabels(), visible=False)
# #plt.savefig('powerlawbin2.png')
# plt.setp(ax2.get_xticklabels(), visible=False)
# plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

# #plot bin 3

# idx_3=np.where((t['RADIUS_PC']>2300)&(t['RADIUS_PC']<3200))
# idxmass_3=np.where((t['MASS_GCORR']>3e5)&(t['RADIUS_PC']>2300)&(t['RADIUS_PC']<3200))

# massdisk3_subset=t['MASS_GCORR'][idxmass_3].data
# fit_3_subset=powerlaw.Fit(massdisk3_subset, xmin=3e5)

# massdisk3=t['MASS_GCORR'][idx_3].data
# fit_3=powerlaw.Fit(massdisk3)


# R3,p3=fit_3_subset.distribution_compare('power_law', 'truncated_power_law')

# ax3=fig.add_subplot(234, sharex=ax1, sharey=ax1)
# fit_3_subset.truncated_power_law.plot_ccdf(label='Trunc. Power Law')
# fit_3_subset.power_law.plot_ccdf(label='Power Law')
# fit_3_subset.plot_ccdf(drawstyle='steps', label='Data')
# plt.axis([10e4, 10e7, 10e-4, 10e-1])
# plt.text(2*10e4, 2*10e-4, '(c)')
# #plt.savefig('powerlawbin3.png')

# #plot bin 4

# idx_4=np.where((t['RADIUS_PC']>3200)&(t['RADIUS_PC']<3900))
# idxmass_4=np.where((t['MASS_GCORR']>3e5)&(t['RADIUS_PC']>3200)&(t['RADIUS_PC']<3900))

# massdisk4_subset=t['MASS_GCORR'][idxmass_4].data
# fit_4_subset=powerlaw.Fit(massdisk4_subset, xmin=3e5)

# massdisk4=t['MASS_GCORR'][idx_4].data
# fit_4=powerlaw.Fit(massdisk4)


# R4,p4=fit_4_subset.distribution_compare('power_law', 'truncated_power_law')

# ax4=fig.add_subplot(235, sharex=ax1, sharey=ax1)
# fit_4_subset.truncated_power_law.plot_ccdf(label='Trunc. Power Law')
# fit_4_subset.power_law.plot_ccdf(label='Power Law')
# fit_4_subset.plot_ccdf(drawstyle='steps', label='Data')
# plt.axis([10e4, 10e7, 10e-4, 10e-1])
# plt.text(2*10e4, 2*10e-4, '(d)')
# plt.setp(ax4.get_yticklabels(), visible=False)
# #plt.savefig('powerlawbin4.png')

# #plot bin 5

# idx_5=np.where((t['RADIUS_PC']>3900)&(t['RADIUS_PC']<4500))
# idxmass_5=np.where((t['MASS_GCORR']>3e5)&(t['RADIUS_PC']>3900)&(t['RADIUS_PC']<4500))

# massdisk5_subset=t['MASS_GCORR'][idxmass_5].data
# fit_5_subset=powerlaw.Fit(massdisk5_subset, xmin=3e5)

# massdisk5=t['MASS_GCORR'][idx_5].data
# fit_5=powerlaw.Fit(massdisk5)

# R5,p5=fit_5_subset.distribution_compare('power_law', 'truncated_power_law')

# ax5=fig.add_subplot(236, sharex=ax1, sharey=ax1)
# fit_5_subset.truncated_power_law.plot_ccdf(label='Trunc. Power Law')
# fit_5_subset.power_law.plot_ccdf(label='Power Law')
# fit_5_subset.plot_ccdf(drawstyle='steps', label='Data')
# plt.axis([10e4, 10e7, 10e-4, 10e-1])
# plt.text(2*10e4, 2*10e-4, '(e)')
# plt.setp(ax5.get_yticklabels(), visible=False)
# #plt.savefig('powerlawbin5.png')
# plt.tight_layout()
# plt.savefig('powerlawsubplots.pdf')

# t=PrettyTable(["Bins", "R", "p", "index", "Truncation mass M_c", "Adamo's index"])
# t.add_row(["0<r<4.5",R,p,-myfit_subset.alpha,1/myfit_subset.truncated_power_law.parameter2 , "N/A"])
# t.add_row(["r<1", -27.555338, 1.139e-13, -1.447356, 19565837, "N/A"])
# t.add_row(["r>1", -46.751947,0.0,-1.731958, 3398336, "N/A"])
# t.add_row(["0<r<0.45", R1, p1, -fit_1_subset.alpha, 1/fit_1_subset.truncated_power_law.parameter2, "N/A"])
# t.add_row(["0.45<r<2.3", R2, p2, -fit_2_subset.alpha,1/fit_2_subset.truncated_power_law.parameter2, -1.90])
# t.add_row(["2.3<r<3.2", R3, p3, -fit_3_subset.alpha,1/fit_3_subset.truncated_power_law.parameter2, -2.20])
# t.add_row(["3.2<r<3.9", R4, p4, -fit_4_subset.alpha,1/fit_4_subset.truncated_power_law.parameter2, -2.20])
# t.add_row(["3.9<r<4.5", R5, p5, -fit_5_subset.alpha,1/fit_5_subset.truncated_power_law.parameter2, -2.70])
# print t

# data=t.get_string()
# with open('powerlawbininfo.txt', 'wb') as f:
#     f.write(data)


