import astropy
import astropy.table
import numpy as np
from astropy.table import Table, Column
import matplotlib.pyplot as plt
from galaxies import Galaxy
import astropy.units as u
import matplotlib as mpl
import astropy.constants as con
mpl.rcParams['font.family']='serif'
mpl.rcParams['font.serif'] = 'Times New Roman'
mygalaxy=Galaxy("M83")

mytable=Table.read('m83.co10.K_props_clfind.fits')
#note that matplotlib understands LaTeX math in the code m

rgal=mygalaxy.radius(ra=(mytable['XPOS']), dec=(mytable['YPOS']))
colrgal=Column(name='RGAL_PC',data=(rgal))
mytable.add_column(colrgal)

condition = np.ones(len(mytable),dtype=np.bool)
condition = (mytable['MAXVAL']/mytable['NOISE']>5)# & (mytable['RADRMS_EXTRAP_DECONV']>20)
idx = np.where(condition)
mytable = mytable[idx]

edge_in = np.array([0,450,2300,3200,3900,4500])
edge_out = np.array([450,2300,3200,3900,4500,10000])

for ins,outs in zip(edge_in,edge_out):
    subset = (mytable['MASS_EXTRAP']>3e5)*(mytable['RGAL_PC']<=outs)*(mytable['RGAL_PC']>ins)
    tt = mytable[subset]
    mass = tt['MASS_EXTRAP'].data
    rho = 3/(4*np.pi)*tt['MASS_EXTRAP']*u.M_sun/\
        (tt['RADRMS_EXTRAP_DECONV']*u.pc)**3
    Pint = rho*(tt['VRMS_EXTRAP_DECONV']*u.km/u.s)**2/con.k_B
    Surfdens = tt['MASS_EXTRAP']*u.M_sun/(tt['RADRMS_EXTRAP_DECONV']*u.pc)**2/np.pi

    tff = ((3*np.pi/(32*con.G*rho))**0.5).to(u.Myr)
    print 'in={0} out={1} Pressure={2:3.2g}, surfdens={3:3.2g}, tff={4:3.2g}'.format(ins,outs,np.nanmedian(Pint.to(u.K/(u.cm**3))),np.nanmedian(Surfdens),np.nanmedian(tff))
