from spectral_cube import SpectralCube, BooleanArrayMask
from astropy.io import fits
import astropy.units as u
import astropy.constants as con
from scipy.special import kv,iv
import numpy as np
from galaxies import Galaxy
from astropy.utils.console import ProgressBar
import scipy.ndimage as nd
import skimage.morphology as morph
from astropy import wcs
from scipy.stats import binned_statistic
from radio_beam import Beam

def lundgren_surfdens(radius):
    # We assume that M83 is 4.8 Mpc away whereas Lundren et al. 2004
    # prefer 4.5 Mpc.  This changes radial scales but not surface
    # densities.
    drat = 4.8/4.5
    # Lundgren et al. use 2.3e20 for Xco, whereas we prefer Xco=2e20
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

def alma_mask():
    pass

def setprogress(axes):
    return(ProgressBar(np.prod(axes)))


def things_vrot(momentname = '/srv/astro/erosolo/m83/data/NGC_5236_RO_MOM1:I:HI:wbb2008.fits'):
    m83 = Galaxy('M83')
    mom1 = fits.getdata(momentname).squeeze()
    hdr = fits.getheader(momentname) 
    x,y = m83.radius(header = hdr,returnXY=True)
    phi = np.arctan2(y.value,x.value)
    r = (x**2+y**2)**0.5
    
    vrot = (mom1/1e3*u.km/u.s-m83.vsys)/(np.sin(m83.inclination)*np.cos(phi))
    idx = np.abs(np.cos(phi))>0.5
    return(r[idx],vrot[idx])
    
def things_profile(momentname = '/srv/astro/erosolo/m83/data/NGC_5236_RO_MOM0:I:HI:wbb2008.fits'):
    m83 = Galaxy('M83')
    mom0 = fits.getdata(momentname).squeeze()
    hdr = fits.getheader(momentname) 
    radius = m83.radius(header = hdr)
    bm = Beam.from_fits_header(hdr)
    nurest = hdr['RESTFREQ']*u.Hz
    radmat = radius.to(u.kpc).value
    dr = 1.0
    nbins=100
    inneredge, outeredge = np.linspace(0,6,nbins), np.linspace(0+dr,6+dr,nbins)
    sdprof = np.zeros(nbins)
    radprof = np.zeros(nbins)
    for ctr, edges in enumerate(zip(inneredge,outeredge)):
        x0,x1 = edges
        idx = (radmat>=x0)*(radmat<x1)
#        import pdb; pdb.set_trace()
        sdprof[ctr] = np.nansum(mom0[idx])/np.sum(np.isfinite(radmat[idx]))
        radprof[ctr] = np.nanmean(radmat[idx])
    sdprof *= np.cos(m83.inclination)*1e-3*bm.jtok(nurest)*0.019
    return radprof,sdprof
        
def alma_profile(momentname = 'm83.moment0.fits'):
    m83 = Galaxy('M83')
    mom0 = fits.getdata(momentname)
    radius = m83.radius(header = fits.getheader(momentname))
    radmat = radius.to(u.kpc).value
    dr = 1.0
    nbins=100
    inneredge, outeredge = np.linspace(0,6,nbins), np.linspace(0+dr,6+dr,nbins)
    sdprof = np.zeros(nbins)
    radprof = np.zeros(nbins)
    for ctr, edges in enumerate(zip(inneredge,outeredge)):
        x0,x1 = edges
        idx = (radmat>=x0)*(radmat<x1)
        sdprof[ctr] = np.nansum(mom0[idx])/np.sum(np.isfinite(radmat[idx]))
        radprof[ctr] = np.nanmean(radmat[idx])
    sdprof *= np.cos(m83.inclination)*4.35
    return radprof, sdprof
        
def alma_moment0(maskname = '/srv/astro/erosolo/m83/data/m83.co10.K_mask.fits',
                 cubename = '/srv/astro/erosolo/m83/data/m83.co10.K_correct.fits'):

    # cube = fits.getdata(cubename)
    # header = fits.getheader(cubename)
    # mask = fits.getdata(maskname)
    # disk = morph.disk(5)
    # disk.shape = (1,)+disk.shape
    # spec = np.ones(3)
    # spec.shape+=(1,1,)
    # mask = nd.binary_dilation(mask,disk*spec)
    # mask *= cube>0
    # cube*=mask
    # moment0 = cube.sum(axis=0)*header['CDELT3']

    #method 2
    sc = SpectralCube.read(cubename, allow_huge_operations=True)
    mask = fits.getdata(maskname)
    disk = morph.disk(5)
    disk.shape = (1,)+disk.shape
    spec = np.ones(3)
    spec.shape+=(1,1,)
    mask = nd.binary_dilation(mask,disk*spec)
    mask = BooleanArrayMask(mask,wcs.WCS(cubename))
    sc = sc.with_mask(mask&(sc>0*u.K))
    moment0 = sc.moment(0)
    

    
def alma_mask():
    emap = fits.getdata('m83.errormap.fits')
    sm_emap = nd.median_filter(emap,footprint=morph.disk(17))
    
def alma_errormap(cube = '/srv/astro/erosolo/m83/data/m83.co10.K_correct.fits'):
    s = SpectralCube.read(cube)
    m83 = Galaxy('M83')
    _,dec,ra = s.world[0,:,:]
    rgal = m83.radius(ra=ra.value,dec=dec.value)
    s2 = s.with_mask(s<0*u.K)
    emap = np.zeros(s.shape[1:])
    pb = setprogress(s.shape[1])
    for x1,x2,slc in s2._iter_rays(0):
        sp = s2.flattened(slc).value
        spdata = sp[np.isfinite(sp)]
        if len(spdata)>0:
#            import pdb; pdb.set_trace()
            emap[x1,x2] = np.median(np.abs(spdata))
        if x2==0:
            pb.update()
    emap*=1.4826 # convert to rms scales
    hdu = fits.PrimaryHDU(emap,header=s.header)
    hdu.writeto('m83.errormap.fits',clobber=True)
