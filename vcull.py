from astropy.io import fits
from astropy.table import Table
import astropy.wcs as wcs
import numpy as np
from scipy.ndimage.interpolation import map_coordinates
from scipy.ndimage.filters import generic_filter
t = Table.read('../measurements/m83.co10.K_props_clfind.fits')

hdu = fits.open('../data/NGC_5236_RO_MOM1:I:HI:wbb2008.fits')
w = wcs.WCS(hdu[0].header)
mom1 = np.squeeze(hdu[0].data)
mom1 = generic_filter(mom1,np.nanmedian,size=(11,11))
x,y,s,v = w.wcs_world2pix(t['XPOS'].data,t['YPOS'].data,np.zeros(len(t)),np.zeros(len(t)),0)

vvals = mom1[y.astype(np.int),x.astype(np.int)]
dv = (vvals/1e3-t['VPOS'].data+22.3)
