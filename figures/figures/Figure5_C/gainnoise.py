import os, sys
lib_path = os.path.abspath(os.path.join('..', 'modules'))
sys.path.append(lib_path)
import fitlorentzian as fl
import getdata as gd
from graph import graph
import numpy as np
import pickle
import matplotlib

dfsw = 'day1_180948_fsw_trace'
fswfreqs, fswpowers, fswparams = gd.getFSWdata(dfsw)
dvna = 'day1_180731_vna_trace'
vnafreqs, vnapowers, vnaphases, vnaparams = gd.getVNAdata(dvna)

cf_res = pickle.load(open('resonance.dat','rb'))
bg = np.abs(cf_res['Renormalisation'])
gainbg_cf = 10* np.log10(bg**2)
centerfreq = cf_res['f0_star']

def amp_cf(freqs, res):
    phases = res['phi0'] - 2*np.arctan(2 * res['Q'] * (freqs/ res['f0_star'] -1))
    z = res['zc'] + res['r0'] * np.exp(1j*phases)
    bg = np.abs(res['Renormalisation'])
    return 10*np.log10(np.abs(z)**2 / bg**2)

noisepar = (fl.fitting(fswfreqs, 10**(fswpowers/10.0), peak=True))
noisebase = noisepar[0]

normgain = vnapowers - gainbg_cf
normnoise = fswpowers - 10*np.log10(noisebase)
midfreq = fswfreqs[len(fswfreqs)/2]

def normfreqs(freqs):
    return 1e6*( freqs - centerfreq)

# fill between
freqscom = np.intersect1d(fswfreqs, vnafreqs)
mask_fsw = np.in1d(fswfreqs, freqscom)
mask_vna = np.in1d(vnafreqs, freqscom)
fbetween = {'x':normfreqs(freqscom), 'y1':normnoise[mask_fsw], 
    'y0':normgain[mask_vna], 'col':'#dfc5fe'}

def inset(fig, ax):
    insax = fig.add_axes([0.15, 0.55, 0.25, 0.35])
    insax.plot(normfreqs(vnafreqs), normgain, color='#ca6641')
    insax.plot(normfreqs(fswfreqs), normnoise, color='#1e488f')
    insax.plot(normfreqs(fswfreqs), amp_cf(fswfreqs, cf_res), 
            color='#8f1402')
    insax.set_xlim([-2,2])
    insax.set_ylim([25,45])
    insax.locator_params(nbins=3)
    insax.fill_between(fbetween['x'], fbetween['y0'], fbetween['y1'], 
            facecolor=fbetween['col'], 
            alpha=0.5, interpolate=True, linewidth=0.0,)

graph(map(normfreqs,[vnafreqs, fswfreqs, fswfreqs]), [normgain, normnoise, amp_cf(fswfreqs, cf_res)],
    'gainnoise.pdf',
    xlab='Frequency offset from $\omega_\mathrm{c}$ (kHz)', 
    ylab='Relative gain and noise (dB)', 
    legends=['$\mathcal{G}$', '$\mathcal{N}$'],
    colors=['#ca6641', '#1e488f', '#8f1402'],
    fbetween=fbetween,
    otherfunc=inset,
    xrange=[-100.,100.],
    yrange=[-3, 45],
    small=True, 
    rast=False)

