from astropy.io import fits
from astropy.table import Table
from galaxies import Galaxy
import matplotlib.pyplot as mpl
mpl.rcParams['font.family']='serif'
mpl.rcParams['font.serif'] = 'Times New Roman'
import aplpy
mygalaxy=Galaxy("M83")

hdu = fits.open('m83.co10.tmax.fits')
rgal=mygalaxy.radius(header=hdu[0].header)/1e3
hdr = hdu[0].header
del hdr['CDELT3']
del hdr['CRPIX3']
del hdr['CRVAL3']
del hdr['NAXIS3']
del hdr['CTYPE3']
del hdr['CUNIT3']
hdr['WCSAXES']=2
newhdu = fits.PrimaryHDU(hdu[0].data,hdr)
radii = fits.PrimaryHDU(rgal,hdr)
map = aplpy.FITSFigure(newhdu)
map.show_colorscale(cmap='Greys',vmin=2.4,vmax=8,stretch='sqrt')
map.add_colorbar()
map.tick_labels.set_font(family='serif')
map.colorbar.set_axis_label_text(r'$T_{\mathrm{max}}\ (\mathrm{K})$')
map.colorbar.set_axis_label_font(family='serif')
map.show_contour(radii,levels=[0.45,2.3,3.2,3.9,4.5],colors=['blue']*5,alpha=0.5,linewidths=[3]*5)

t = Table.read('m83.co10.K_props_cprops.fits')
map.show_markers(t['XPOS'],t['YPOS'],edgecolors=['yellow']*len(t),facecolors=['black']*len(t),s=4)
map.save('m83_tmax.pdf')
