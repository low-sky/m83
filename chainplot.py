import numpy as np
from astropy.table import Table
from astropy.table import Column
import matplotlib.pyplot as plt
from galaxies import Galaxy
import astropy.units as u
import matplotlib as mpl
import plfit_mcstat as plf
import scipy.stats as ss
import seaborn as sns

mpl.rcParams['font.family']='serif'
mpl.rcParams['font.serif'] = 'Times New Roman'
mpl.rcParams['font.size']=12

tin = Table.read('m83.co10.K_props_clfind.fits')
mygalaxy = Galaxy("M83")
minmass_global = 3.35e5
minmass = 3.35e5
# find cloud's galactocentric distance
rgal = mygalaxy.radius(ra=(tin['XPOS']), dec=(tin['YPOS']))
# fit1,fit2 = 'truncated_power_law','schechter'
# fit1,fit2 = 'power_law','truncated_power_law'
# fit1,fit2 = 'power_law','schechter'
# #add those distances to the fits table
colrgal = Column(name='RGAL_PC', data=(rgal))
tin.add_column(colrgal)

edge_in = np.array([3200])
edge_out = np.array([3900])
type = 'trunc'
for ins, outs in zip(edge_in, edge_out):
    subset = (tin['MASS_GCORR'] > 0) * (tin['RGAL_PC'] <= outs) *\
        (tin['RGAL_PC'] > ins)
    tt = tin[subset]
    mass = np.sort(tt['MASS_GCORR'].data)
    minmass = np.max([minmass_global])
    mass = mass[mass > minmass_global]
    optresult = plf.plfit_ksstat(mass / minmass, type=type, method='ad')
    sampler = plf.plfit_emcee(mass / minmass, type=type, method='ad',
                              steps=5000, thin=5)
    pvec = optresult.x
    g = sns.jointplot(sampler.flatchain[:,0] - 1,
                      sampler.flatchain[:,1] + np.log10(3.35e5),
                      stat_func=None,
                      kind="hex").set_axis_labels(r"$\beta$",
                                                  r"$\log_{10}(M_{c}/M_{\odot})$")
    g.fig.set_size_inches(4.0,4.25)
    g.ax_joint.set_xlim(-2.5,1.0)
    g.ax_joint.set_ylim(5.5,8.5)
plt.tight_layout()
plt.savefig('samplerdist.pdf')
plt.close()
plt.clf()
