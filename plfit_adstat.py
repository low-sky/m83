import numpy as np
import scipy.optimize as opt
import emcee
from mpmath import gammainc, mpf
from numpy import vectorize
gammainc = vectorize(gammainc)
mpfv = vectorize(mpf)


def cdf_truncpl(x, alpha=-1.0, xmin=1, xmax=np.inf):
    xvals = np.sort(x) / xmax
    xvals = mpfv(xvals.tolist())
    cdf_theoretical = 1 - gammainc(alpha, a=xvals) /\
        gammainc(alpha, a=(xmin / xmax))
    return(cdf_theoretical.astype('float'))


def adstat_trunc(pvec, data):
    alpha = pvec[0]
    xmax = 1e1**pvec[1]
#    xmin = 1.0
    sorted_xvals = np.sort(data)
    npts = np.float(len(sorted_xvals))
    cdf_theoretical = np.zeros_like(sorted_xvals)
#    for idx, elt in enumerate(sorted_xvals):
#        cdf_theoretical[idx] = 1 - gammainc(alpha, elt / xmax) /\
#            gammainc(alpha, (xmin / xmax))
    cdf_theoretical = cdf_truncpl(sorted_xvals, alpha=alpha, xmax=xmax)
    cdf_empirical = 1 - np.linspace(1 / npts / 2,
                                    1 - 1 / npts / 2, npts)
    adstat = len(sorted_xvals) * \
        np.sum((cdf_empirical - cdf_theoretical)**2 /
               (cdf_theoretical * (1 - cdf_theoretical)))
    print (alpha, xmax, adstat)
    return(adstat)


def adstat_bounded(pvec, data):
    alpha = pvec[0]
    xmax = 1e1**pvec[1]
    xmin = 1.0
    sorted_xvals = np.sort(data)
    npts = np.float(len(sorted_xvals))
    if np.any(data > xmax):
        return(np.inf)
    cdf_theoretical = 1 - (1 - (xmin / sorted_xvals)**alpha) /\
        (1 - (xmin / xmax)**alpha)
    cdf_empirical = 1 - np.linspace(1 / npts / 2,
                                    1 - 1 / npts / 2, npts)
#    forwards = np.log(cdf_theoretical)
#    backwards = np.log(1 - cdf_theoretical)
#    adstatlog = -npts - (1 / npts) * np.sum((forwards + backwards) *
#                                            (2 * np.linspace(1,
#                                                             npts, npts) - 1))
    adstat = len(sorted_xvals) * \
        np.sum((cdf_empirical - cdf_theoretical)**2 /
               (cdf_theoretical * (1 - cdf_theoretical)))
#    print (pvec, adstat)
    return(adstat)


def pareto_rng(n=500, xmin=1.0, xmax=10.0, alpha=1.0):
    y = np.random.rand(n)
#    x = xmin / (1 - y)**(1 / alpha)
    xbounded = ((1 - y * (1 - (xmin / xmax)**alpha)) /
                xmin**alpha)**(-1 / alpha)
    return(xbounded)


def plfit_adstat(indata, xmin=1.0):
    data = indata / xmin  # Normalize to range
    result = opt.minimize(adstat_trunc,
                          np.array([-1.0, np.log10(data.max()) + 0.1]),
                          bounds=[(-1.5, -0.3),
                                  (np.log10(data.max()),
                                   np.log10(data.max()) + 2)],
                          args=(data),
                          method='Nelder-Mead')

    return(result)


def logprob_plfit(pvec, data):
    if (pvec[0] < 0.0) or (pvec[0] > 1.5):
        return(-np.inf)
    if (pvec[1] > 20):
        return(-np.inf)
    adstat = adstat_trunc(pvec, data)
    logprob = -4 * ((np.log(adstat / data.size) + 0.25)**2 / 2 / 0.67**2)
    return logprob


def grid_eval(indata, xmin=1.0):
    opt_result = plfit_adstat(indata, xmin=xmin)
    if opt_result.success:
        alpha, logxmax = np.meshgrid(opt_result.x[0] *
                                     np.linspace(0.5, 1.5, 101),
                                     opt_result.x[1] *
                                     np.linspace(0.95, 2.05, 101))

        adstatmat = np.zeros(alpha.size)
        for idx, p in enumerate(zip(alpha.ravel(),
                                    logxmax.ravel())):
            adstatmat[idx] = adstat_bounded(p, indata)
        adstatmat.shape = alpha.shape
        return(adstatmat, alpha, logxmax)


def plfit_emcee(indata, xmin=1.0):
    moddata  = indata / xmin
    moddata = moddata[moddata >= 1.0]
    opt_result = plfit_adstat(moddata, xmin=1.0)
    print(opt_result.x)
    if opt_result.success:
        ndim = 2
        nwalkers = ndim * 10
        p0 = np.zeros((nwalkers, ndim))
        p0[:, 0] = np.random.randn(nwalkers) * 0.1 + opt_result.x[0]
        p0[:, 1] = np.abs(np.random.randn(nwalkers) * 0.01) + opt_result.x[1]
        sampler = emcee.EnsembleSampler(nwalkers,
                                        ndim, logprob_plfit,
                                        args=[moddata])
        pos, prob, state = sampler.run_mcmc(p0, 200)
        sampler.reset()
        pos, prob, state = sampler.run_mcmc(pos, 2000)
#        import pdb; pdb.set_trace()
        return sampler
