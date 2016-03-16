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
import scipy.interpolate as interp


def channelShift(x,ChanShift):
    # Shift a spectrum by a set number of channels.  
    ftx = np.fft.fft(x)
    m = np.fft.fftfreq(len(x))
    phase = np.exp(2*np.pi*m*1j*ChanShift)
    x2 = np.real(np.fft.ifft(ftx*phase))
    return(x2)
    
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

# for M83 V_helio - 10.919938 = V_LSRK
def cloud_veldisp(catalog ='/home/erosolow/bigdata/erosolo/m83/measurements/m83.co10_props_cprops.fits'):
    m83 = Galaxy('M83')
    radius = Galaxy.radius(ra=t['RA'],dec=t['DEC'])

def ism_hidisp(momentname = '/home/erosolow/bigdata/erosolo/m83/data/NGC_5236_RO_MOM1:I:HI:wbb2008.fits',
               cubename = '/home/erosolow/bigdata/erosolo/m83/data/NGC_5236_RO_CUBE:I:HI:wbb2008.fits.gz',
               dr = 0.25, nbins=100):
    s = SpectralCube.read(cubename)
    mom1 = fits.getdata(momentname)
    m83 = Galaxy('M83')
    radius = m83.radius(header=s.header)
    intfunc = interp.interp1d(s.spectral_axis.value,np.arange(s.shape[0]),bounds_error = False)
    channel = intfunc(mom1.squeeze().ravel())
    channel.shape = s.shape[1:]
    inneredge, outeredge = np.linspace(0,6,nbins), np.linspace(0+dr,6+dr,nbins)
    vaxis = s.spectral_axis.to(u.km/u.s).value
    sigma = np.zeros(nbins)
    for ctr, edges in enumerate(zip(inneredge,outeredge)):
        x0,x1 = edges
        idx = (radius>=x0*u.kpc)*(radius<x1*u.kpc)
        ymat, xmat = np.indices(s.shape[1:])
        dchan = channel[idx]-s.shape[0]/2
        accumspec = np.zeros(s.shape[0])

        for deltachan,yspec,xspec in zip(dchan,ymat[idx],xmat[idx]):
            if np.isfinite(deltachan):
                accumspec += channelShift(s[:,yspec,xspec],deltachan)
        import pdb; pdb.set_trace()
        labels,_ = nd.measurements.label(accumspec>0)
        pk = np.argmax(accumspec)
        roi = (labels == labels[pk])
        v0 = np.sum(accumspec[roi]*vaxis[roi])/np.sum(accumspec[roi])
        sigma[ctr] = (np.sum(accumspec[roi]*(vaxis[roi]-v0)**2)/np.sum(accumspec[roi]))**0.5
    return sigma

    
def ism_veldisp(momentname = '/home/erosolow/bigdata/erosolo/m83/data/m83.mom1.fits',
                cubename = '/home/erosolow/bigdata/erosolo/m83/data/m83.co10.K_correct.fits',
                dr = 0.25, nbins=100):
    sc = SpectralCube.read(cubename)
    _, dec, ra = scworld[0,:,:]
    mom1 = fits.getdata(momentname)
#    mom1 = nd.generic_filter(mom1,np.nanmedian,size=1)
    wcs_mom1 = wcs.WCS(momentname)
    xvals, yvals = wcs_mom1.celestial.wcs_world2pix(ra,dec,0)
    vvals = nd.interpolation.map_coordinates(mom1.squeeze(),[yvals.ravel(),xvals.ravel()],order=1)*u.m/u.s
    intfunc = interp.interp1d(sc.spectral_axis.value,np.arange(scshape[0]),bounds_error = False)
    channel = intfunc(vvals.to(u.km/u.s).value-10.919938)
    channel.shape = sc.shape[1:]
    m83 = Galaxy('M83')
    radius = m83.radius(header=sc.header)
    inneredge, outeredge = np.linspace(0,6,nbins), np.linspace(0+dr,6+dr,nbins)
    vaxis = sc.spectral_axis.to(u.km/u.s).value
    sigma = np.zeros(nbins)
    for ctr, edges in enumerate(zip(inneredge,outeredge)):
        x0,x1 = edges
        idx = (radius>=x0*u.kpc)*(radius<x1*u.kpc)
        ymat, xmat = np.indices(sc.shape[1:])
        dchan = channel[idx]-50
        accumspec = np.zeros(sc.shape[0])
        for deltachan,yspec,xspec in zip(dchan,ymat[idx],xmat[idx]):
            if np.isfinite(deltachan):
                accumspec += channelShift(sc[:,yspec,xspec],deltachan)
        
        labels,_ = nd.measurements.label(accumspec>0)
        pk = np.argmax(accumspec)
        roi = (labels == labels[pk])
        v0 = np.sum(accumspec[roi]*vaxis[roi])/np.sum(accumspec[roi])
        sigma[ctr] = (np.sum(accumspec[roi]*(vaxis[roi]-v0)**2)/np.sum(accumspec[roi]))**0.5
    return sigma
    
def things_vrot(momentname = '/home/erosolow/bigdata/erosolo/m83/data/m83.mom1.fits',
                dr = 0.25,nbins=100):
    m83 = Galaxy('M83')
    mom1 = fits.getdata(momentname).squeeze()
    mom1 = nd.median_filter(mom1,size=1)
    hdr = fits.getheader(momentname) 
    x,y = m83.radius(header = hdr,returnXY=True)
    phi = np.arctan2(y.value,x.value)
    r = (x**2+y**2)**0.5
    r = r.to(u.kpc).value
    #Optical frame Vsys
    vrot = (mom1*u.km/u.s-506*u.km/u.s)/(np.sin(m83.inclination)*np.cos(phi))
    idx = np.abs(np.cos(phi))>0.25
    inneredge, outeredge = np.linspace(0,6,nbins), np.linspace(0+dr,6+dr,nbins)
    vprof = np.zeros(nbins)
    radprof = np.zeros(nbins)
    vscatter = np.zeros(nbins)
    vrot = vrot.to(u.km/u.s).value
    for ctr, edges in enumerate(zip(inneredge,outeredge)):
        x0,x1 = edges
        idx = (r>=x0)*(r<x1)
        vprof[ctr] = np.nanmedian(vrot[idx])
        vscatter[ctr] = np.nanmedian(np.abs(vrot[idx]-vprof[ctr]))*1.4826
        radprof[ctr] = np.nanmean(r[idx])
    return radprof,vprof,vscatter
    
def things_profile(momentname = '/home/erosolow/bigdata/erosolo/m83/data/NGC_5236_RO_MOM0:I:HI:wbb2008.fits',
                   dr=0.5,nbins=100):
    m83 = Galaxy('M83')
    mom0 = fits.getdata(momentname).squeeze()
    hdr = fits.getheader(momentname) 
    radius = m83.radius(header = hdr)
    bm = Beam.from_fits_header(hdr)
    nurest = hdr['RESTFREQ']*u.Hz
    radmat = radius.to(u.kpc).value
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
        
def alma_profile(momentname = 'm83.moment0.fits',
                 dr=0.5,nbins=100):
    m83 = Galaxy('M83')
    mom0 = fits.getdata(momentname)
    radius = m83.radius(header = fits.getheader(momentname))
    radmat = radius.to(u.kpc).value
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
        
def alma_moment0(maskname = '/home/erosolow/bigdata/erosolo/m83/data/m83.co10.K_mask.fits',
                 cubename = '/home/erosolow/bigdata/erosolo/m83/data/m83.co10.K_correct.fits'):

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
    moment1 = sc.moment(1)
    return moment0,moment1

def alma_mask():
    emap = fits.getdata('m83.errormap.fits')
    sm_emap = nd.median_filter(emap,footprint=morph.disk(17))
    
def alma_errormap(cube = '/home/erosolow/bigdata/erosolo/m83/data/m83.co10.K_correct.fits'):
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
