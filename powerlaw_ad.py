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
import plfit_adstat as plf

minmass = 4e5

#  get info about m83
mygalaxy = Galaxy("M83")
print(mygalaxy)

#  load fits file
t = Table.read('m83.co10.K_props_clfind.fits')

#  find cloud's galactocentric distance
rgal = mygalaxy.radius(ra=(t['XPOS']), dec=(t['YPOS']))

#  add those distances to the fits table
colrgal = Column(name='RADIUS_PC', data=(rgal))
t.add_column(colrgal)
#  print(t)

#  create a loop to input multiple bin noundaries
inneredge = np.array([0, 450, 2300, 3200, 3900])
outeredge = np.array([450, 2300, 3200, 3900, 4500])

column_names = ['Inner edge (pc)', 'Outer edge (pc)',
                'GMC index', 'R', 'p',
                'Truncation mass (M$_\odot$)',
                'Largest cloud (M$_\odot$)',
                '5th largest cloud (M$_\odot$)']
column_types = ['f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4']
table = Table(names=column_names, dtype=column_types)

for inneredge, outeredge in zip(inneredge, outeredge):
    idx = np.where((t['RADIUS_PC'] >= inneredge) &
                 (t['RADIUS_PC'] < outeredge))
    mass = t['MASS_GCORR'][idx].data

    optresult = plf.plfit_adstat(mass / minmass)
    pvec = optresult.x
#  don't have to create an index for the xmin mass - defined in the fit_subset
    table.add_row()
    table[-1]['Inner edge (pc)'] = inneredge
    table[-1]['Outer edge (pc)'] = outeredge
    table[-1]['GMC index'] = pvec[0]
    table[-1]['R'] = 0
    table[-1]['p'] = 0
    table[-1]['Truncation mass (M$_\odot$)'] = pvec[1] * minmass
    table[-1]['Largest cloud (M$_\odot$)'] = mass.max()
    table[-1]['5th largest cloud (M$_\odot$)'] =\
        np.sort(t['MASS_GCORR'][idx])[-5]

    print(table)

table.write('m83bininfo_ad.fits', overwrite=True)
