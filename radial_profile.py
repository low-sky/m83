from spectral_cube import SpectralCube
from astropy.io import fits
import astropy.units as u
import astropy.constants as con
from scipy.special import kv,iv
import numpy as np
def lundgren_surfdens(radius):
    # We assume that M83 is 4.8 Mpc away whereas Lundren et al. 2004
    # prefer 4.5 Mpc.  This changes radial scales but not surface
    # densities.
    drat = 4.8/4.5


    # Lundren et al. use 2.3e20 for Xco, whereas we prefer Xco=2e20
    xrat = 2.0/2.3 
    sigma_h2 = xrat * (216*np.exp(-(radius/(drat*649*u.pc))**2)+
                       120*np.exp(-(radius/(drat*2265*u.pc))))*u.M_sun/u.pc**2
    return sigma_h2

def lundgren_vrot(radius):
    diskmass = 6e10*u.M_sun*(4.8/4.5)
    diskradius = (2.7*u.kpc*4.8/4.5)
    y = (radius/diskradius).to(u.dimensionless_unscaled)
    y = y.value
    vc = (2*con.G*diskmass*y**2/(diskradius)*
          (kv(0,y)*iv(0,y)-kv(1,y)*iv(1,y)))**0.5
    vc = vc.to(u.km/u.s)
    return vc

