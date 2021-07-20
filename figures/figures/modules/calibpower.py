import numpy as np
import fitlorentzian as fl
from math import pi

# Parameters of the system used for the noise calibration (see SI VIII)
alpha = -1.6 # attenuation
# HEMT level
HEMT_bg = -128.59 # in dBm/Hz, from average of fitted noise curves
HEMT_T = 4.27 # in Kelvin
kB = 1.38e-23; hbar = 1.05457e-34
HEMT_n = kB * HEMT_T / (hbar * 2 * np.pi * 4.125e9) # in quanta/sec

def FSWpower_to_HEMTin(FSWpow, subtractHEMT=True, HEMT_bg=HEMT_bg, 
        HEMT_n=HEMT_n, **keywargs):
    """
    FSWpow in dBm/Hz

    returns equivalent power at input with HEMT noise subtracted and in photons/second
    """
    Pin = HEMT_n * 10**((FSWpow - HEMT_bg)/10.0)
    if subtractHEMT:
        Pin = Pin - HEMT_n
    return  Pin

def FSWpower_to_Saa(FSWpow, subtractHEMT=True, HEMT_bg=HEMT_bg, alpha=alpha, **keywargs):
    """
    FSWpow in dBm/Hz

    returns equivalent power at ouput of device with HEMT noise subtracted and in photons/second
    """
    Pin = FSWpower_to_HEMTin(FSWpow, subtractHEMT=subtractHEMT, HEMT_bg=HEMT_bg, **keywargs)
    Saa = Pin / 10.0**(alpha/10.0)
    return Saa

def FSWpower_to_DUTin(freqs, FSWpow, gainf, alpha=alpha, subtractHEMT=True, **keywargs):
    """
    FSWpow: in dBm/Hz
    gainf: function of gain (power, lin.) as function of frequency
    alpha: attenuation (if negative) in dB

    returns equivalent power at input with HEMT noise subtracted and in photons/second
    """
    gaincurve = gainf(freqs)
    FSWpow_out = FSWpower_to_HEMTin(FSWpow, subtractHEMT=subtractHEMT, **keywargs) # returns linear power
    FSWpow_out = FSWpow_out / gaincurve
    # attenuation increases noise at input of device
    FSWpow_out = FSWpow_out / 10.0**(alpha/10.0)
    return FSWpow_out

def FSWpower_to_neff(freqs, FSWpow, gainf, gain, coop, eta, alpha=alpha, ncav=0.0, **keywargs):
    """
    FSWpow: in dBm/Hz, assumes no probe (max=max of noise)
    gainf: function of gain (power, lin.) as function of frequency
    gain: gain on resonnance (to subtract HEMT noise)
    coop: cooperativity (such that kappa_eff = (1-coop) kappa
    eta: kappa_ex / kappa
    alpha: attenuation (if negative) in dB

    returns mechanical occupations (equivalent thermal) to account for given input noise in resonance
    """
    FSWpow_in = FSWpower_to_DUTin(freqs, FSWpow, gainf, alpha=alpha, subtractHEMT=False, **keywargs) # returns linear power
    # remove probe tone
    indprobe = np.argmax(FSWpow_in)
    npars = (fl.fitting(np.delete(freqs, indprobe), 
        np.delete(FSWpow_in, indprobe),
        peak=False))
    N_omc = npars[0] + (npars[1]/pi) / (npars[3]/2.0)
    # subract HEMT noise
    N_omc = N_omc - ( HEMT_n / 10**(gain/10.0) )/ 10.0**(alpha/10.0)
    neff = (N_omc * (1-2*eta-coop)**2 - 4*eta*(1-eta)*(ncav+0.5)) / (4*eta*coop) - 0.5
    return neff, N_omc, npars

def Nadded(coop, eta, neff):
    return (4*coop*eta*(neff+0.5) + 4*eta*(1-eta)*0.5) / (1-coop-2*eta)**2

if __name__ == "__main__":
    print 'standard quantum limit = {} phot/ s Hz'.format(Nadded(1,eta,0))
    print 'N added by HEMT = {} phot/ s Hz'.format(HEMT_n)

