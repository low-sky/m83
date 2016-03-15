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

#get info about m83
mygalaxy=Galaxy("M83")
print(mygalaxy)

#load fits file
t=Table.read('/Users/erik/m83.co10.K_props_cprops.fits')

#find cloud's galactocentric distance
rgal=mygalaxy.radius(ra=(t['XPOS']), dec=(t['YPOS']))

#add those distances to the fits table
colrgal=Column(name='RADIUS_PC',data=(rgal))
t.add_column(colrgal)
#print(t)

#create a loop to input multiple bin noundaries
inneredge=np.array([0, 450, 2300, 3200, 3900])
outeredge=np.array([450, 2300, 3200, 3900, 4500])

for inneredge, outeredge in zip(inneredge, outeredge):
    idx=np.where((t['RADIUS_PC']>=inneredge)&(t['RADIUS_PC']<outeredge))
    mass=t['MASS_EXTRAP'][idx].data
    #don't have to create an index for the xmin mass - defined in the fit_subset
    fit=powerlaw.Fit(mass)
    fit_subset=powerlaw.Fit(mass, xmin=3e5)
    R,p=fit_subset.distribution_compare('power_law', 'truncated_power_law')
    R2,p2 = fit_subset.distribution_compare('power_law', 'schechter')
    import pdb; pdb.set_trace()
    print(-fit.alpha, -fit_subset.alpha, R, p, 1/fit_subset.truncated_power_law.parameter2)   

#table=PrettyTable(["alpha", "subset alpha", "p", "index", "Truncation mass M_c"])
#table.add_row([-fit.alpha, -fit_subset.alpha,R,p,1/myfit_subset.truncated_power_law.parameter2])
    
#data=table.get_string()
#with open('powerlawloopbininfo.txt', 'wb') as f:
 #   f.write(data)
        
    
