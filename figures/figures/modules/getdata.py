import os.path as path
import numpy as np
import glob

def bname(d):
    return path.basename(d)

def get_setdat(dirname, setdat='dat', warning=True):
    ls = glob.glob(dirname + '/*.' + setdat)
    if warning:
        if len(ls)>1:
            print('More than one {} file!').format(setdat)
    return ls[0]

def get_a_b(dirname, instrumentname, propertyname, datatype=float):
    basename = path.basename(dirname)
    setfile = get_setdat(dirname, 'set')
    with open(setfile, 'r') as setfile:
        for line in setfile:
            if 'Instrument: {}'.format(instrumentname) in line:
                for line_inst in setfile:
                    if propertyname in line_inst:
                        prop = datatype(line_inst.split()[1])
                        return prop

def getsourceP(dirname, srcst):
    power = get_a_b(dirname, srcst, 'power')
    return power

def getsourceF(dirname, srcst):
    frequency = get_a_b(dirname, srcst, 'frequency')
    return frequency

def getfswBW(dirname, fswname='fsw'):
    bd = get_a_b(dirname, fswname, 'bandwidth')
    return bd


def getFSWdata(d, srcst='src1', dBmHz=True, warning=True):
    """
    Returns data from FSW sweeps:
        freqs: list of frequencies in GHz
        fswpowers: list of powers in dBm/Hz
        params: dictionary with parameters (rbw, pumppower)
    """
    datfile = get_setdat(d, 'dat', warning=warning)
    data = np.loadtxt(datfile)

    bw = getfswBW(d)
    p = getsourceP(d, srcst)

    freqs = data[:,0]

    fswpowers = data[:,1]
    # noise spectrum in dBm / Hz
    if dBmHz:
        fswpowers = fswpowers - 10*np.log10(bw*1.0645)

    params = {
            'rbw':bw,
            'pumppower':p,
            }
    
    return freqs, fswpowers, params

def getVNAdata(d, srcst='src1', warning=True):
    """
    Returns data from VNA sweeps:
        freqs: list of frequencies in GHz
        vnapowers: list of powers in dB
        phases: list of phases in radians
        params: dictionary with parameters (rbw, pumppower)
    """
    datfile = get_setdat(d, 'dat', warning=warning)
    data = np.loadtxt(datfile)

    p = getsourceP(d, srcst)

    freqs = data[:,0]
    vnapowers = data[:,1]
    phases = data[:,3]

    params = {
            'pumppower':p,
            }
    
    return freqs, vnapowers, phases, params

def get_S_from_FSW(d, probesrc='src3'):
    """
    Returns Signal peak from FSW data:
        d: FSW data directory
        probesrc: src for the probe signal
    returns
        fprobe: frequency in GHz
        Spow: gain in arbitrary units (dB)
    """
    datfile = get_setdat(d, 'dat')
    data = np.loadtxt(datfile)

    freqs = data[:,0]
    fswpowers = data[:,1]

    fprobe = getsourceF(d, probesrc)
    fprobe = fprobe / 1e9

    Spow = max( fswpowers[np.logical_and( fprobe-1e-6 < freqs, freqs < fprobe+1e-6 )] )
    return fprobe, Spow

def getgaindat(swpd):
    gaindatf = swpd + '/gaindata.dat'
    if path.exists(gaindatf):
        data = np.loadtxt(gaindatf)
        freqs = data[:,0]
        Spow = data[:,1]
    else:
        dirsfsw = sorted(glob(swpd + '/*fsw_trace'))
        n = len(dirsfsw)
        freqs = np.zeros(n-1)
        Spow = np.zeros(n-1)
        for i in range(n-1): # drop first element
            f, p = gd.get_S_from_FSW(dirsfsw[i+1], 'src3')
            f = f
            freqs[i] = f
            Spow[i] = p
            print(f, p)

        # reorder data
        sortedind = np.argsort(freqs)
        freqs = freqs[sortedind]
        Spow = Spow[sortedind]
        comment = 'Gain data of {}\nfreqs (GHz), gain (dB, arbitrary)'.format(gd.bname(swpd))
        np.savetxt(gaindatf, np.transpose((freqs, Spow)),
                header=comment)
    return freqs, Spow

