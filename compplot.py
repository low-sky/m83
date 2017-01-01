import matplotlib.pyplot as plt
from astropy.table import Table
import glob
import numpy as np
import matplotlib as mpl

mpl.rcParams['font.family']='serif'
mpl.rcParams['font.serif'] = 'Times New Roman'
mpl.rcParams['font.size']=12

filelist = glob.glob('/Users/erik/code/m83/comprun*fits')

cfrac = np.zeros((51, len(filelist)))
for i, thisfile in enumerate(filelist):
    t = Table.read(thisfile)
    cfrac[:, i] = t['FracRecover'] / t['FracRecover'].max()
    lm = t['logmass']
compfrac = 1 / (1 + np.exp(-14.8 * (lm - 5.525)))

mnvec = np.mean(cfrac, axis=1)
sdvec = np.std(cfrac, axis=1)
fig, ax = plt.subplots(1, 1)
fig.set_size_inches(4.25, 3.5)
ax.semilogx(1e1**lm, np.mean(cfrac, axis=1), label='Completeness')
ax.plot(1e1**lm, compfrac, color='red', linestyle='--', label='Model')
ax.fill_between(1e1**lm, mnvec - sdvec,
                mnvec + sdvec, alpha=0.25)
ax.set_xlabel(r'$M_{\mathrm{GMC}}\ (M_{\odot})$')
ax.set_ylabel('Recovery Fraction')
ax.set_ylim(0, 1)
ax.hlines(0.5, 1e5, 1e6, lw=2, alpha=0.25, color='grey')
ax.vlines(3.35e5, 0, 1, lw=2, alpha=0.25, color='grey')
ax.legend(loc='upper left', fontsize=10)
fig.tight_layout()
plt.savefig('compfrac.pdf')
plt.close()
plt.clf()
