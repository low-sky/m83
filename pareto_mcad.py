import numpy as np


def adstat_bounded(data, alpha=1.0, xmin=1.0, xmax=100):
    if np.any(data > xmax):
        return(np.inf)
    sorted_xvals = np.sort(data)
    npts = np.float(len(sorted_xvals))
    cdf_theoretical = 1 - (1 - (xmin / sorted_xvals)**alpha) /\
        (1 - (xmin / xmax)**alpha)
    cdf_empirical = 1 - np.linspace(0.5 / npts, 1 - 0.5 / npts, npts)
    adstat = len(sorted_xvals) * np.sum((cdf_empirical -
                                         cdf_theoretical)**2 /
                                        (cdf_theoretical *
                                         (1 - cdf_theoretical)))
    forwards = np.log(cdf_theoretical)
    backwards = np.log(1 - cdf_theoretical[::-1])
    adstatlog = -npts - (1 / npts) * np.sum((forwards + backwards) *
                                            (2 * np.linspace(1,
                                                             npts,
                                                             npts) - 1))
    return(adstat, adstatlog)


npts = np.logspace(1.5, 3.5, 15)
xmaxarr = np.logspace(1, 4, 15)
alphaarr = np.linspace(0.5, 1.5, 10)

ntrials = 100
xmin = 1.0

adstat = np.zeros((npts.size, xmaxarr.size, alphaarr.size, ntrials))
adstatlog = np.zeros((npts.size, xmaxarr.size, alphaarr.size, ntrials))

for i1, n in enumerate(npts):
    for i2, xmax in enumerate(xmaxarr):
        for i3, aaa in enumerate(alphaarr):
            for t in np.arange(ntrials):

                y = np.random.rand(n)
                x = xmin / (1 - y)**(1 / aaa)
                xbounded = ((1 - y * (1 - (xmin / xmax)**aaa)) /
                            xmin**aaa)**(-1 / aaa)
                adstat[i1, i2, i3, t], tmp = adstat_bounded(xbounded,
                                                            alpha=aaa,
                                                            xmax=xmax)
                adstatlog[i1, i2, i3, t] = tmp
