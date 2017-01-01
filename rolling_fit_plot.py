# import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table

t = Table.read('rolling_plfit.fits')

fig, ax  = plt.subplots(1,1)
fig.set_size_inches(4.25,3.5)
plt.semilogy(t['Rmean'], t['Mtrun_trunc'], 'ro')
plt.fill_between(t['Rmean'], t['Mtrun15_trunc'],
                 t['Mtrun85_trunc'], 'ro', color='grey')
plt.savefig('rollingfit.pdf')
