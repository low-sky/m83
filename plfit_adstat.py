import numpy as np
import scipy.optimize as opt

def adstat_bounded(pvec, data):
    alpha = pvec[0]
    xmax = pvec[1]
    xmin = 1.0
    sorted_xvals = np.sort(data)
    npts = np.float(len(sorted_xvals))
    if np.any(data > xmax):
        return(np.inf)
    cdf_theoretical = 1-(1 - (xmin/sorted_xvals)**alpha)/(1-(xmin/xmax)**alpha)
    cdf_empirical = 1 - np.linspace(1/npts/2, 1-1/npts/2, npts)
    forwards = np.log(cdf_theoretical)
    backwards = np.log(1-cdf_theoretical)
    adstatlog = -npts - (1/npts) * np.sum((forwards+backwards)*
                                          (2*np.linspace(1,npts,npts)-1))
    adstat = len(sorted_xvals)*np.sum((cdf_empirical-cdf_theoretical)**2/
                                      (cdf_theoretical*(1-cdf_theoretical)))

    print (pvec, adstatlog)
    return(adstat)


def plfit_adstat(indata, xmin = 1.0):
    data = indata /  xmin  # Normalize to range
    result = opt.minimize(adstat_bounded, 
                          np.array([1.0, 4*data.max()]), 
                          bounds = [(0.5,1.5),(0.5*data.max(),10*data.max())],
                          args = (data))
    return(result)
