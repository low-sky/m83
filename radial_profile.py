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
from astropy.table import Table, Column
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['font.family']='serif'

from scipy.signal import savgol_filter


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
    y = (0.5*radius/diskradius).to(u.dimensionless_unscaled)
    y = y.value
    vc = (2*con.G*diskmass*y**2/(diskradius)*
          (kv(0,y)*iv(0,y)-kv(1,y)*iv(1,y)))**0.5
    vc = vc.to(u.km/u.s)
    return vc

def alma_mask():
    pass

def setprogress(axes):
    return(ProgressBar(np.prod(axes)))

def critical_masses():
    t = Table.read('m83.profiles.fits')
    sd = (t['Surfdens_HI'].data*t['Surfdens_HI'].unit+\
         t['Surfdens_H2'].data*t['Surfdens_H2'].unit)

    veldisp = ((t['sigma_v_hi'].data*t['sigma_v_hi'].unit)**2+\
              (t['sigma_v_h2'].data*t['sigma_v_h2'].unit)**2)**0.5


#    veldisp = t['sigma_v_hi'].data*t['sigma_v_hi'].unit
#    idx = t['sigma_v_h2']>t['sigma_v_hi'] 
#    veldisp[idx] = t['sigma_v_h2'][idx].data*t['sigma_v_h2'][idx].unit
    t['sigma_v_h2'].data*t['sigma_v_h2'].unit
    vrot = t['Vrot'].data*t['Vrot'].unit
    radius = t['Radius'].data*t['Radius'].unit

    Mjeans = (veldisp**4/(np.pi*con.G**2*sd)).to(u.M_sun)
    Ltoomre = (np.pi*con.G*sd/(vrot/radius)**2).to(u.pc)
    Mtoomre = (sd*np.pi*Ltoomre**2).to(u.M_sun)

    r4d =np.insert(radius,0,0*u.kpc)
    
    v4d = np.insert(vrot,0,0*u.km/u.s)
    grad = np.gradient(v4d*r4d)
    grad = grad[1:]
    import pdb; pdb.set_trace()
    kappa = (2*vrot/radius**2*grad/
             (radius[1]-radius[0]))**0.5


    diskmass = 6e10*u.M_sun*(4.8/4.5)
    diskradius = (2.7*u.kpc*4.8/4.5)
    y = (0.5*radius/diskradius).to(u.dimensionless_unscaled)
    y=y.value
    kappa2 = 4*con.G*diskmass/diskradius**3*\
             ((2*iv(0,y)+y*iv(1,y))*kv(0,y)-
              (y*iv(0,y)+iv(1,y))*kv(1,y))

    Ltoomre_full = 4*np.pi**2*con.G*sd/kappa2#**2
    Mtoomre_full = (sd*np.pi*Ltoomre_full**2).to(u.M_sun)

    inneredge = np.array([0,0.45,2.3,3.2,3.9])*u.kpc
    outeredge = np.array([0.45,2.3,3.2,3.9,4.5])*u.kpc

    for (xin,xout) in zip(inneredge,outeredge):
        idx = (radius >= xin)*(radius<xout)
        mj = np.mean(Mjeans[idx])
        mt = np.mean(Mtoomre[idx])
        mtf = np.mean(Mtoomre_full[idx])
        print "{0:3.1f} & {1:3.1f} & {2:3.1f} \\".format(mj.value/1e6,mt.value/1e6,mtf.value/1e6)

    import pdb; pdb.set_trace()


def allthethings():
    radius_co, surfdens_co = alma_profile()
    radius_hi, surfdens_hi = things_profile()
    radius_rc, vrot, vscatter = things_vrot()
    veldmol = ism_veldisp()
    c0 = Column(radius_co*u.kpc,name='Radius')
    c1 = Column(surfdens_co*u.M_sun/u.pc**2,name='Surfdens_H2')
    c2 = Column(surfdens_hi*u.M_sun/u.pc**2,name='Surfdens_HI')
    c3 = Column(vrot*u.km/u.s,name='Vrot')
    c4 = Column(vscatter*u.km/u.s,name='V_scatter')
    c5 = Column(veldmol*u.km/u.s,name='sigma_v_hi')
    t = Table()
    for c in [c0,c1,c2,c3,c4,c5]:
        t.add_column(c)
    t.write('m83.profiles.fits',overwrite=True)

    veldism = ism_hidisp()
    c6 = Column(veldism*u.km/u.s,name='sigma_v_h2')
    for c in [c6]:
        t.add_column(c)
    t.write('m83.profiles.fits',overwrite=True)




    import pdb; pdb.set_trace()

# for M83 V_helio - 10.919938 = V_LSRK
def cloud_veldisp(catalog ='/home/erosolow/fastdata/m83/measurements/m83.co10_props_cprops.fits'):
    m83 = Galaxy('M83')
    radius = Galaxy.radius(ra=t['RA'],dec=t['DEC'])

def ism_hidisp(momentname = '/home/erosolow/fastdata/m83/data/NGC_5236_RO_MOM1:I:HI:wbb2008.fits',
               cubename = '/home/erosolow/fastdata/m83/data/NGC_5236_RO_CUBE:I:HI:wbb2008.fits.gz',
               dr = 0.25, nbins=100):
    sc = SpectralCube.read(cubename,allow_huge_operations=True,
                           mode='readonly',memmap=False)
    sc = sc[4:-4,:,:]
    mom1 = fits.getdata(momentname)
    m83 = Galaxy('M83')
    radius = m83.radius(header=sc.header)
    intfunc = interp.interp1d(sc.spectral_axis.value,np.arange(sc.shape[0]),bounds_error = False)
    channel = intfunc(mom1.squeeze().ravel())
    channel.shape = sc.shape[1:]
    inneredge, outeredge = np.linspace(0,6,nbins), np.linspace(0+dr,6+dr,nbins)
    vaxis = sc.spectral_axis.to(u.km/u.s).value
    sigma = np.zeros(nbins)
    for ctr, edges in enumerate(zip(inneredge,outeredge)):
        x0,x1 = edges
        idx = (radius>=x0*u.kpc)*(radius<x1*u.kpc)
        ymat, xmat = np.indices(sc.shape[1:])
        dchan = channel[idx]-sc.shape[0]/2
        accumspec = np.zeros(sc.shape[0])

        for deltachan,yspec,xspec in zip(dchan,ymat[idx],xmat[idx]):
            if np.isfinite(deltachan):
                accumspec += channelShift(sc[:,yspec,xspec],deltachan)
        labels,_ = nd.measurements.label(accumspec>0)
        pk = np.argmax(accumspec)
        roi = (labels == labels[pk])
        v0 = np.sum(accumspec[roi]*vaxis[roi])/np.sum(accumspec[roi])
        sigma[ctr] = (np.sum(accumspec[roi]*(vaxis[roi]-v0)**2)/np.sum(accumspec[roi]))**0.5
        print(ctr)
    return sigma

    
def ism_veldisp(momentname = '/home/erosolow/fastdata/m83/data/m83.mom1.fits',
                cubename = '/home/erosolow/fastdata/m83/data/m83.co10.K_correct.fits',
                dr = 0.25, nbins=100):
    sc = SpectralCube.read(cubename,allow_huge_operations=True,
                           mode='readonly',memmap=False)
    _, dec, ra = sc.world[0,:,:]
    mom1 = fits.getdata(momentname)
#    mom1 = nd.generic_filter(mom1,np.nanmedian,size=1)
    wcs_mom1 = wcs.WCS(momentname)
    xvals, yvals = wcs_mom1.celestial.wcs_world2pix(ra,dec,0)
    vvals = nd.interpolation.map_coordinates(mom1.squeeze(),[yvals.ravel(),xvals.ravel()],order=1)*u.m/u.s
    intfunc = interp.interp1d(sc.spectral_axis.value,np.arange(sc.shape[0]),bounds_error = False)
    channel = intfunc(vvals.to(u.km/u.s).value-10.919938)
    channel.shape = sc.shape[1:]
    m83 = Galaxy('M83')
    radius = m83.radius(header=sc.header)
    inneredge, outeredge = np.linspace(0,6,nbins), np.linspace(0+dr,6+dr,nbins)
    vaxis = sc.spectral_axis.to(u.km/u.s).value
    sigma = np.zeros(nbins)
    for ctr, edges in enumerate(zip(inneredge[0:10],outeredge[0:10])):
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
        print ctr, sigma[ctr]
    return sigma
    
def things_vrot(momentname = '/home/erosolow/fastdata/m83/data/NGC_5236_RO_MOM1:I:HI:wbb2008.fits',
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
    vrot = (mom1/1e3*u.km/u.s-506*u.km/u.s)/(np.sin(m83.inclination)*np.cos(phi))
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
    
def things_profile(momentname = '/home/erosolow/fastdata/m83/data/NGC_5236_RO_MOM0:I:HI:wbb2008.fits',
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
        
def alma_moment0(maskname = '/home/erosolow/fastdata/m83/data/m83.co10.K_mask.fits',
                 cubename = '/home/erosolow/fastdata/m83/data/m83.co10.K_correct.fits'):

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
    
def alma_errormap(cube = '/home/erosolow/fastdata/m83/data/m83.co10.K_correct.fits'):
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

def profile_plots():
    import matplotlib.pyplot as plt

    t = Table.read('m83.profiles.fits')

    sd = t['Surfdens_HI'].data*t['Surfdens_HI'].unit+\
         t['Surfdens_H2'].data*t['Surfdens_H2'].unit
    veldisp = ((t['sigma_v_hi'].data*t['sigma_v_hi'].unit)**2+\
              (t['sigma_v_h2'].data*t['sigma_v_h2'].unit)**2)**0.5

    rad = t['Radius'].data
    vrot = t['Vrot'].data
    plt.clf()



    plt.figure(figsize=(4.5,6))
    ax0 = plt.subplot(311)
    plt.semilogy(rad, t['Surfdens_H2'].data,label=r'H$_2$',linestyle=':')
    plt.plot(rad, t['Surfdens_HI'].data,label=r'HI',linestyle='--')
    plt.plot(rad, t['Surfdens_H2'].data+t['Surfdens_HI'].data,label='Total',linestyle='-')
    plt.plot(rad, lundgren_surfdens(rad*u.kpc).value,label = 'L04',alpha=0.7,lw=4)
    plt.setp(ax0.get_xticklabels(), visible=False)
    plt.legend(loc=1,fontsize=9)
    plt.ylabel(r'$\Sigma\ (M_{\odot}\ \mathrm{pc}^{-2})$')



    ax1 = plt.subplot(312,sharex=ax0)
    plt.plot(rad,vrot,label='THINGS')
    plt.plot(rad,(lundgren_vrot(rad*u.kpc).to(u.km/u.s)).value,
             label='L04',linestyle='--')
    plt.fill_between(rad,vrot+t['V_scatter'].data/2,vrot-t['V_scatter']/2,color='0.75')
    plt.legend(loc=4,fontsize=10)
    plt.ylabel(r'V (km s$^{-1}$)')
    plt.setp(ax1.get_xticklabels(), visible=False)
    
    ax2 = plt.subplot(313,sharex=ax0)
    plt.plot(rad,savgol_filter(veldisp.value,11,3,mode='nearest'))
    plt.ylim([0,40])
    plt.xlabel(r'R (kpc)')
    plt.ylabel(r'$\sigma_v$ (km s$^{-1}$)')



    plt.tight_layout()
    plt.savefig('profiles1.pdf')
    plt.clf()
    plt.figure(figsize=(5,3.5))

    plt.savefig('surfdens.pdf')




    pass
