# import matplotlib.pyplot as plt
from astropy.table import Table
import glob

filelist = glob.glob('~/Desktop/complim/comprun*fits')

for thisfile in filelist:
    t = Table.read(thisfile)

