import numpy as np
import scipy.optimize as opt
import emcee
from mpmath import gammainc, mpf
from numpy import vectorize
import scipy.stats as ss
import scipy.special as scispec
gammainc = vectorize(gammainc)
mpfv = vectorize(mpf)


def gammaincplus2(s, x):
    vals = scispec.gamma(s + 2) * scispec.gammaincc(s + 2, x)
    vals = (vals - x**(s + 1) *
            np.exp(-x) - (s + 1) * x**s * np.exp(-x))
    vals /= (s * (s + 1))
    return(vals)


def cdf_truncpl(x, alpha, xmin, xmax):  # alpha=-1.0, xmin=1, xmax=np.inf):
    xvals = np.sort(x / xmin) / xmax
    if alpha <= -2:
        xvals = mpfv(xvals.tolist())
        cdf_theoretical = 1 - gammainc(alpha, a=xvals) /\
            gammainc(alpha, a=(xmin / xmax))
        return(cdf_theoretical.astype('float'))
    else:
        cdf_theoretical = 1 - gammaincplus2(alpha, xvals) /\
            gammaincplus2(alpha, 1.0 / xmax)
        return(cdf_theoretical)


def cdf_boundpl(x, alpha, xmin, xmax):
    sorted_xvals = np.sort(x / xmin)
    cdf_theoretical = (1 - (xmin / sorted_xvals)**(-alpha)) /\
        (1 - (xmin / xmax)**(-alpha))
    return cdf_theoretical


def ksstat_trunc(pvec, data, returnlnprob=False):
    alpha = pvec[0]
    xmax = 1e1**pvec[1]
    sorted_xvals = np.sort(data)
    D, p = ss.kstest(sorted_xvals, cdf_truncpl,
                     args=(alpha, 1, xmax), mode='asymp')
    if returnlnprob:
        if np.isnan(np.log(p)):
            return(-np.inf)
        return 2 * np.log(p)
    return(D)


def ksstat_bounded(pvec, data, returnlnprob=False):
    alpha = pvec[0]
    xmax = 1e1**pvec[1]
    sorted_xvals = np.sort(data)
    D, p = ss.kstest(sorted_xvals, cdf_boundpl,
                     args=(alpha, 1, xmax), mode='asymp')
    if returnlnprob:
        return 2 * np.log(p)
    return(D)


def ksstat_schechter(pvec, data, returnlnprob=False):
    alpha = -0.999  # pvec[0]
    xmax = 1e1**pvec[1]
    sorted_xvals = np.sort(data)
    D, p = ss.kstest(sorted_xvals, cdf_truncpl,
                     args=(alpha, 1, xmax), mode='asymp')
    if returnlnprob:
        if np.isnan(p):
            return(-np.inf)
        return 2 * np.log(p)
    return(D)


def ksstat_purepl(pvec, data, returnlnprob=False):
    alpha = pvec[0]
    xmax = 1e20
    sorted_xvals = np.sort(data)
    D, p = ss.kstest(sorted_xvals, cdf_truncpl,
                     args=(alpha, 1, xmax), mode='asymp')
    if returnlnprob:
        return 2 * np.log(p)
    return(D)


def adstat_bounded(pvec, data, returnlnprob=False):
    alpha = pvec[0]
    xmax = 1e1**pvec[1]
    xmin = 1.0
    sorted_xvals = np.sort(data)
    npts = np.float(len(sorted_xvals))
    cdf_theoretical = cdf_boundpl(sorted_xvals, alpha, xmin, xmax)
    cdf_empirical = np.linspace(1 / npts / 2,
                                1 - 1 / npts / 2, npts)
    adstat = len(sorted_xvals) * \
        np.sum((cdf_empirical - cdf_theoretical)**2 /
               (cdf_theoretical * (1 - cdf_theoretical)))
    if returnlnprob:
        logprob = -data.size * (np.log(adstat / data.size) +
                                0.25)**2 / 2 / 0.67**2
        if adstat < 0:
            return(-np.inf)
        return(-adstat / 2)
    return(adstat)


def adstat_schechter(pvec, data, returnlnprob=False):
    alpha = -0.999  # pvec[0]
    xmax = 1e1**pvec[1]
    xmin = 1.0
    sorted_xvals = np.sort(data)
    npts = np.float(len(sorted_xvals))
    cdf_theoretical = cdf_truncpl(sorted_xvals, alpha, xmin, xmax)
    cdf_empirical = np.linspace(1 / npts / 2,
                                1 - 1 / npts / 2, npts)
    adstat = len(sorted_xvals) * \
        np.sum((cdf_empirical - cdf_theoretical)**2 /
               (cdf_theoretical * (1 - cdf_theoretical)))
    if returnlnprob:
        if adstat < 0:
            return(-np.inf)
        return(-adstat / 2)
    return(adstat)


def adstat_trunc(pvec, data, returnlnprob=False):
    alpha = pvec[0]
    xmax = 1e1**pvec[1]
    xmin = 1.0
    sorted_xvals = np.sort(data)
    npts = np.float(len(sorted_xvals))
    cdf_theoretical = cdf_truncpl(sorted_xvals, alpha, xmin, xmax)
    cdf_empirical = np.linspace(1 / npts / 2,
                                1 - 1 / npts / 2, npts)
    adstat = len(sorted_xvals) * \
        np.sum((cdf_empirical - cdf_theoretical)**2 /
               (cdf_theoretical * (1 - cdf_theoretical)))
    if returnlnprob:
        if np.isnan(adstat):
            return(-np.inf)
        return(-adstat / 2)
    return(adstat)


def pareto_rng(n=500, xmin=1.0, xmax=10.0, alpha=1.0):
    y = np.random.rand(n)
    xbounded = ((1 - y * (1 - (xmin / xmax)**alpha)) /
                xmin**alpha)**(-1 / alpha)
    return(xbounded)


def plfit_ksstat(indata, xmin=1.0, type='bound', method='ks'):
    data = indata / xmin  # Normalize to range
    if method == 'ks':
        if type == 'bound':
            fitstat = ksstat_bounded
        if type == 'trunc':
            fitstat = ksstat_trunc
        if type == 'sch':
            fitstat = ksstat_schechter
        if type == 'purepl':
            fitstat = ksstat_purepl

    if method == 'ad':
        if type == 'bound':
            fitstat = adstat_bounded
        if type == 'sch':
            fitstat = adstat_schechter
        if type == 'trunc':
            fitstat = adstat_trunc

    result = opt.minimize(fitstat,
                          np.array([-0.7, np.log10(data.max()) + 0.1]),
                          bounds=[(-1.5, -0.3),
                                  (np.log10(data.max()),
                                   np.log10(data.max()) + 2)],
                          args=(data),
                          method='Nelder-Mead')

    return(result)


def logprob_plfit(pvec, data, type='bound', method='ks'):
    if (pvec[1] > 20):
        return(-np.inf)
    if method == 'ks':
        if type == 'bound':
            logprob = ksstat_bounded(pvec, data, returnlnprob=True)
        if type == 'trunc':
            logprob = ksstat_trunc(pvec, data, returnlnprob=True)
        if type == 'sch':
            logprob = ksstat_schechter(pvec, data, returnlnprob=True)
        if type == 'purepl':
            logprob = ksstat_purepl(pvec, data, returnlnprob=True)

    if method == 'ad':
        if type == 'sch':
            logprob = adstat_schechter(pvec, data, returnlnprob=True)
        if type == 'bound':
            logprob = adstat_bounded(pvec, data, returnlnprob=True)
        if type == 'trunc':
            logprob = adstat_trunc(pvec, data, returnlnprob=True)
    return logprob


def grid_eval(indata, xmin=1.0):
    opt_result = plfit_ksstat(indata, xmin=xmin)
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


def plfit_emcee(indata, xmin=1.0, type='bound', method='ks'):
    moddata  = indata / xmin
    opt_result = plfit_ksstat(moddata, xmin=1.0,
                              type=type, method=method)
    print(opt_result.x)
    if opt_result.success:
        ndim = 2
        nwalkers = ndim * 10
        p0 = np.zeros((nwalkers, ndim))
        p0[:, 0] = np.random.randn(nwalkers) * 0.1 + opt_result.x[0]
        p0[:, 1] = np.random.randn(nwalkers) * 0.1 + opt_result.x[1]
        sampler = emcee.EnsembleSampler(nwalkers,
                                        ndim, logprob_plfit,
                                        args=[moddata],
                                        kwargs={'type': type,
                                                'method': method})
        pos, prob, state = sampler.run_mcmc(p0, 200)
        sampler.reset()
        pos, prob, state = sampler.run_mcmc(pos, 2000, thin=10)

        return sampler
