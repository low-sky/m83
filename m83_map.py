from astropy.io import fits
import aplpy
from galaxies import Galaxy
import matplotlib as mpl
mpl.rcParams['font.family']='serif'
mpl.rcParams['font.serif'] = 'Times'
mygalaxy=Galaxy("M83")

hdu = fits.open('m83.co10.tmax.fits')
rgal=mygalaxy.radius(header=hdu[0].header)
hdr = hdu[0].header
del hdr['CDELT3']
del hdr['CRPIX3']
del hdr['CRVAL3']
del hdr['NAXIS3']
del hdr['CTYPE3']
del hdr['CUNIT3']
hdr['WCSAXES']=2
newhdu = fits.PrimaryHDU(hdu[0].data,hdr)
map = aplpy.FITSFigure(newhdu)
map.show_colorscale(cmap='Greys',vmin=2,vmax=10)
map.add_colorbar()
map.tick_labels.set_font(family='serif',size='small')
map.colorbar.set_axis_label_text(r'$T_{\mathrm{max}}\ (\mathrm{K})$')
map.colorbar.set_axis_label_font(family='serif',size='small')
map.save('m83_tmax.pdf')
