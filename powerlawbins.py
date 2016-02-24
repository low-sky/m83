# -*- coding: utf-8 -*-

#import libraries
import astropy
import astropy.table
import numpy as np
from astropy.table import Table
from astropy.table import Column
import matplotlib.pyplot as plt
import powerlaw
from galaxies import Galaxy
import astropy.units as u
from tabulate import tabulate
from prettytable import PrettyTable

#get info about m83
mygalaxy=Galaxy("M83")
print(mygalaxy)

#load fits file
t=Table.read('/home/pafreema/Documents/m83.co10_props_cprops.fits')

#find cloud's galactocentric distance
rgal=mygalaxy.radius(ra=(t['XPOS']), dec=(t['YPOS']))

#add those distances to the fits table
colrgal=Column(name='RADIUS_PC',data=(rgal))
t.add_column(colrgal)
#print(t)

#fit for the whole galaxy
mass=t['MASS_EXTRAP'].data
myfit=powerlaw.Fit(mass)

#plot bin 1

#index the file so you get clouds within a certain radius
idx_1=np.where(t['RADIUS_PC']<450)
idxmass_1=np.where((t['MASS_EXTRAP']>myfit.xmin)&(t['RADIUS_PC']<450))

#pull out the mass variable and fit it
massdisk1=t['MASS_EXTRAP'][idx_1].data
fit_1=powerlaw.Fit(massdisk1)
massdisk1_subset=t['MASS_EXTRAP'][idxmass_1].data
fit_1_subset=powerlaw.Fit(massdisk1_subset, xmin=myfit.xmin)



#fit_1.plot_ccdf() #look at the data with the fit
print(fit_1.alpha) #returns the fitted index
print(fit_1_subset.alpha)

#compares two different forms of the distribution
#R<0 means the data is more consistent with the right distribution, p is the significance of the result
R1,p1=fit_1_subset.distribution_compare('power_law', 'truncated_power_law')
print(R1,p1)

plt.figure(1)
fit_1_subset.truncated_power_law.plot_ccdf(label='Trunc. Power Law')
fit_1_subset.power_law.plot_ccdf(label='Power Law')
fit_1_subset.plot_ccdf(drawstyle='steps', label='data')
plt.xlabel(r'$Mass\ M_{\odot}$')
plt.ylabel(r'$N\ (>M)$')
plt.legend(loc=3)
plt.savefig('powerlawbin1.png')

#plot bin 2

idx_2=np.where((t['RADIUS_PC']>450)&(t['RADIUS_PC']<2300))
idxmass_2=np.where((t['MASS_EXTRAP']>myfit.xmin)&(t['RADIUS_PC']>450)&(t['RADIUS_PC']<2300))

massdisk2_subset=t['MASS_EXTRAP'][idxmass_2].data
fit_2_subset=powerlaw.Fit(massdisk2_subset, xmin=myfit.xmin)
massdisk2=t['MASS_EXTRAP'][idx_2].data
fit_2=powerlaw.Fit(massdisk2)

#fit_1.plot_ccdf() #look at the data with the fit
print(fit_2.alpha) #returns the fitted index
print(fit_2_subset.alpha)

R2,p2=fit_2_subset.distribution_compare('power_law', 'truncated_power_law')
print(R2,p2)

plt.figure(2)
fit_2_subset.truncated_power_law.plot_ccdf(label='Trunc. Power Law')
fit_2_subset.power_law.plot_ccdf(label='Power Law')
fit_2_subset.plot_ccdf(drawstyle='steps', label='data')
plt.xlabel(r'$Mass\ M_{\odot}$')
plt.ylabel(r'$N\ (>M)$')
plt.legend()
plt.savefig('powerlawbin2.png')

#plot bin 3

idx_3=np.where((t['RADIUS_PC']>2300)&(t['RADIUS_PC']<3200))
idxmass_3=np.where((t['MASS_EXTRAP']>2496477.5)&(t['RADIUS_PC']>2300)&(t['RADIUS_PC']<3200))

massdisk3_subset=t['MASS_EXTRAP'][idxmass_3].data
fit_3_subset=powerlaw.Fit(massdisk3_subset, xmin=2496477.5)

massdisk3=t['MASS_EXTRAP'][idx_3].data
fit_3=powerlaw.Fit(massdisk3)

#fit_1.plot_ccdf() #look at the data with the fit
print(fit_3.alpha) #returns the fitted index
print(fit_3_subset.alpha)

R3,p3=fit_3_subset.distribution_compare('power_law', 'truncated_power_law')
print(R3,p3)

plt.figure(3)
fit_3_subset.truncated_power_law.plot_ccdf(label='Trunc. Power Law')
fit_3_subset.power_law.plot_ccdf(label='Power Law')
fit_3_subset.plot_ccdf(drawstyle='steps', label='data')
plt.xlabel(r'$Mass\ M_{\odot}$')
plt.ylabel(r'$N\ (>M)$')
plt.legend()
plt.savefig('powerlawbin3.png')

#plot bin 4

idx_4=np.where((t['RADIUS_PC']>3200)&(t['RADIUS_PC']<3900))
idxmass_4=np.where((t['MASS_EXTRAP']>myfit.xmin)&(t['RADIUS_PC']>3200)&(t['RADIUS_PC']<3900))

massdisk4_subset=t['MASS_EXTRAP'][idxmass_4].data
fit_4_subset=powerlaw.Fit(massdisk4_subset, xmin=myfit.xmin)

massdisk4=t['MASS_EXTRAP'][idx_4].data
fit_4=powerlaw.Fit(massdisk4)

#fit_1.plot_ccdf() #look at the data with the fit
print(fit_4.alpha) #returns the fitted index
print(fit_4_subset.alpha)

R4,p4=fit_4_subset.distribution_compare('power_law', 'truncated_power_law')
print(R4,p4)

plt.figure(4)
fit_4_subset.truncated_power_law.plot_ccdf(label='Trunc. Power Law')
fit_4_subset.power_law.plot_ccdf(label='Power Law')
fit_4_subset.plot_ccdf(drawstyle='steps', label='data')
plt.xlabel(r'$Mass\ M_{\odot}$')
plt.ylabel(r'$N\ (>M)$')
plt.legend(loc=3)
plt.savefig('powerlawbin4.png')

#plot bin 5

idx_5=np.where((t['RADIUS_PC']>3900)&(t['RADIUS_PC']<4500))
idxmass_5=np.where((t['MASS_EXTRAP']>myfit.xmin)&(t['RADIUS_PC']>3900)&(t['RADIUS_PC']<4500))

massdisk5_subset=t['MASS_EXTRAP'][idxmass_5].data
fit_5_subset=powerlaw.Fit(massdisk5_subset, xmin=myfit.xmin)

massdisk5=t['MASS_EXTRAP'][idx_5].data
fit_5=powerlaw.Fit(massdisk5)

#fit_1.plot_ccdf() #look at the data with the fit
print(fit_5.alpha) #returns the fitted index
print(fit_5_subset.alpha)

R5,p5=fit_5_subset.distribution_compare('power_law', 'truncated_power_law')
print(R5,p5)

plt.figure(5)
fit_5_subset.truncated_power_law.plot_ccdf(label='Trunc. Power Law')
fit_5_subset.power_law.plot_ccdf(label='Power Law')
fit_5_subset.plot_ccdf(drawstyle='steps', label='data')
plt.xlabel(r'$Mass\ M_{\odot}$')
plt.ylabel(r'$N\ (>M)$')
plt.legend(loc=3)
plt.savefig('powerlawbin5.png')

t=PrettyTable(["Bins", "R", "p", "index", "Adamo's index"])
t.add_row(["0<r<4.5", -0.680348, 0.243417, -2.425428, "N/A"])
t.add_row(["r<1", -3.470554, 0.008424, -1.958932, "N/A"])
t.add_row(["r>1", -1.066727, 0.144116, -3.205864, "N/A"])
t.add_row(["0<r<0.45", R1, p1, -fit_1_subset.alpha, "N/A"])
t.add_row(["0.45<r<2.3", R2, p2, -fit_2_subset.alpha, -1.90])
t.add_row(["2.3<r<3.2", R3, p3, -fit_3_subset.alpha, -2.20])
t.add_row(["3.2<r<3.9", R4, p4, -fit_4_subset.alpha, -2.20])
t.add_row(["3.9<r<4.5", R5, p5, -fit_5_subset.alpha, -2.70])
print t
    
data=t.get_string()
with open('powerlawbininfo.txt', 'wb') as f:
    f.write(data)
        
    