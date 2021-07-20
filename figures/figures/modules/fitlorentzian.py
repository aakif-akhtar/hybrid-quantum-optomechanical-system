import numpy as np
from scipy.optimize import leastsq
from scipy.optimize import curve_fit

def lorentzian(x, p):
    # parameters
    # constant, area, center, width
    # assumes that it is given np.Array. If not, problems with +, etc
    return p[0] + (p[1]/ np.pi) * (np.abs(p[3])/2.0) / ((x- p[2])**2+(0.5 * p[3])**2 )

def lorentzian_alt(x, p0, p1, p2, p3):
    # parameters
    # constant, area, center, width
    # assumes that it is given np.Array. If not, problems with +, etc
    p = [p0, p1, p2, p3]
    return lorentzian(x, p)

def amplitude(p):
    return (p[1]/np.pi) / (p[3]/2.0)

def maximum(p):
    return p[0] + (p[1]/np.pi) / (p[3]/2.0)

def residual_func(fitlog=False):
    def residuals(p, y, x):
        if fitlog:
            return np.log(np.abs(y/lorentzian(x,p)))
        else:
            return y - lorentzian(x, p)
    return residuals

def normalize(freq, linpow, lp0, norm):
    f0 = 0.0; fscale=1.0; Pscale=1.0
    if normalize:
        # get data closer to 1 for better precision
        f0 = freq[len(freq)//2]
        lfreq = freq - f0
        fscale = 10**round(np.log10(abs(lfreq[-1]))) 
        lfreq = lfreq/fscale
        Pscale = 10**round(np.log10(abs(linpow[0])))
        llinpow = linpow/Pscale
        lp0[0] = lp0[0]/Pscale
        lp0[1] = lp0[1]/(Pscale * fscale)
        lp0[2] = (lp0[2]-f0)/fscale
        lp0[3] = lp0[3]/fscale
    else:
        lfreq = freq
        llinpow = linpow

    return lfreq, llinpow, lp0, Pscale, fscale, f0

def haspeak(linpow):
    base = 0.5*(linpow[-1] + linpow[0])
    l = len(linpow)
    inside = linpow[l/3: 2*l/3]
    outside = linpow[0:l/3]
    peaktrue = ((np.log(inside)).mean() - np.log(outside).mean()) > 0
    return peaktrue

def fitting(freq, linpow, p0 = [0.1, 1.0, 0.1, 1.0], 
        norm=True, guess=True, peak=None, fitlog=False,
        returnerror=False):
    if fitlog and returnerror:
        print ("Cannot fit logscale and return error at this point!")
    if guess:
        if peak is None:
            peak = haspeak(linpow) 
        lp0 = [0,0,0,0]
        if peak:
            lp0[0] = min(linpow)
            lp0[3] = (freq[-1]-freq[0])/3.  # photography 1/3 rule
            lp0[1] = (max(linpow)-min(linpow)) * lp0[3] * (np.pi / 2.)
            lp0[2] = freq[np.argmax(linpow)]
        else:
            lp0[0] = max(linpow)
            lp0[3] = (freq[-1]-freq[0])/3.  # photography 1/3 rule
            lp0[1] = (max(linpow)-min(linpow)) * lp0[3] * (np.pi / 2.)
            lp0[2] = freq[np.argmin(linpow)]
    else:
        lp0 = list(p0) # deep copy
    lfreq, llinpow, lp0, Pscale, fscale, f0 = normalize(freq, linpow, lp0, norm)

    if not returnerror:
        pbest = leastsq(residual_func(fitlog), lp0, args=(llinpow, lfreq), 
                full_output=1)
        bestpar = pbest[0]
    else:
        bestpar, pcov = curve_fit(lorentzian_alt, lfreq, llinpow, lp0)
        std_dev = np.sqrt(np.diag(pcov))
        std_dev = [
                std_dev[0] * Pscale,
                std_dev[1] * Pscale * fscale,
                std_dev[2] * fscale,
                std_dev[3] * fscale] 
    # rescaling :
    bestpar = [
            bestpar[0] * Pscale,
            bestpar[1] *Pscale * fscale,
            bestpar[2] * fscale + f0,
            np.abs(bestpar[3]) * fscale] 
    if not returnerror:
        return bestpar
    else:
        return bestpar, std_dev

