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

#plot bin 1

#index the file so you get clouds within a certain radius
idx_1=np.where(t['RADIUS_PC']<450)

#pull out the mass variable and fit it
massdisk1=t['MASS_EXTRAP'][idx_1].data
fit_1=powerlaw.Fit(massdisk1)

#fit_1.plot_ccdf() #look at the data with the fit
print(fit_1.alpha) #returns the fitted index

#compares two different forms of the distribution
#R<0 means the data is more consistent with the right distribution, p is the significance of the result
R1,p1=fit_1.distribution_compare('power_law', 'truncated_power_law')
print(R1,p1)

plt.figure(1)
fit_1.truncated_power_law.plot_ccdf(label='Trunc. Power Law')
fit_1.power_law.plot_ccdf(label='Power Law')
fit_1.plot_ccdf(drawstyle='steps', label='data')
plt.xlabel(r'$Mass\ M_{\odot}$')
plt.ylabel(r'$N\ (>M)$')
plt.legend()
plt.savefig('powerlawbin1.png')

#plot bin 2

idx_2=np.where((t['RADIUS_PC']>450)&(t['RADIUS_PC']<2300))

massdisk2=t['MASS_EXTRAP'][idx_2].data
fit_2=powerlaw.Fit(massdisk2)

#fit_1.plot_ccdf() #look at the data with the fit
print(fit_2.alpha) #returns the fitted index

R2,p2=fit_2.distribution_compare('power_law', 'truncated_power_law')
print(R2,p2)

plt.figure(2)
fit_2.truncated_power_law.plot_ccdf(label='Trunc. Power Law')
fit_2.power_law.plot_ccdf(label='Power Law')
fit_2.plot_ccdf(drawstyle='steps', label='data')
plt.xlabel(r'$Mass\ M_{\odot}$')
plt.ylabel(r'$N\ (>M)$')
plt.legend()
plt.savefig('powerlawbin2.png')

#plot bin 3

idx_3=np.where((t['RADIUS_PC']>2300)&(t['RADIUS_PC']<3200))

massdisk3=t['MASS_EXTRAP'][idx_3].data
fit_3=powerlaw.Fit(massdisk3)

#fit_1.plot_ccdf() #look at the data with the fit
print(fit_3.alpha) #returns the fitted index

R3,p3=fit_3.distribution_compare('power_law', 'truncated_power_law')
print(R3,p3)

plt.figure(3)
fit_3.truncated_power_law.plot_ccdf(label='Trunc. Power Law')
fit_3.power_law.plot_ccdf(label='Power Law')
fit_3.plot_ccdf(drawstyle='steps', label='data')
plt.xlabel(r'$Mass\ M_{\odot}$')
plt.ylabel(r'$N\ (>M)$')
plt.legend()
plt.savefig('powerlawbin3.png')

#plot bin 4

idx_4=np.where((t['RADIUS_PC']>3200)&(t['RADIUS_PC']<3900))

massdisk4=t['MASS_EXTRAP'][idx_4].data
fit_4=powerlaw.Fit(massdisk4)

#fit_1.plot_ccdf() #look at the data with the fit
print(fit_4.alpha) #returns the fitted index

R4,p4=fit_4.distribution_compare('power_law', 'truncated_power_law')
print(R4,p4)

plt.figure(4)
fit_4.truncated_power_law.plot_ccdf(label='Trunc. Power Law')
fit_4.power_law.plot_ccdf(label='Power Law')
fit_4.plot_ccdf(drawstyle='steps', label='data')
plt.xlabel(r'$Mass\ M_{\odot}$')
plt.ylabel(r'$N\ (>M)$')
plt.legend()
plt.savefig('powerlawbin4.png')

#plot bin 5

idx_5=np.where((t['RADIUS_PC']>3900)&(t['RADIUS_PC']<4500))

massdisk5=t['MASS_EXTRAP'][idx_5].data
fit_5=powerlaw.Fit(massdisk5)

#fit_1.plot_ccdf() #look at the data with the fit
print(fit_5.alpha) #returns the fitted index

R5,p5=fit_5.distribution_compare('power_law', 'truncated_power_law')
print(R5,p5)

plt.figure(5)
fit_5.truncated_power_law.plot_ccdf(label='Trunc. Power Law')
fit_5.power_law.plot_ccdf(label='Power Law')
fit_5.plot_ccdf(drawstyle='steps', label='data')
plt.xlabel(r'$Mass\ M_{\odot}$')
plt.ylabel(r'$N\ (>M)$')
plt.legend()
plt.savefig('powerlawbin5.png')


t=PrettyTable(["Bins", "R", "p", "index", "Adamo's index"])
t.add_row(["0<r<4.5", -0.610329, 0.269231, 2.436146, "N/A"])
t.horizontal_char
t.add_row(["r<1", -0.406056, 0.367496, 2.541322, "N/A"])
t.add_row(["r>1", -0.644785, 0.256127, 3.309108, "N/A"])
t.horizontal_char
t.add_row(["0<r<0.45", R1, p1, -fit_1.alpha, "N/A"])
t.add_row(["0.45<r<2.3", R2, p2, -fit_2.alpha, -1.90])
t.add_row(["2.3<r<3.2", R3, p3, -fit_3.alpha, -2.20])
t.add_row(["3.2<r<3.9", R4, p4, -fit_4.alpha, -2.20])
t.add_row(["3.9<r<4.5", R5, p5, -fit_5.alpha, -2.70])
print t
    
data=t.get_string()
with open('powerlawbininfo.txt', 'wb') as f:
    f.write(data)
        
    