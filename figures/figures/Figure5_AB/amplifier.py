import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pickle
import scipy.optimize as scipopt
import os, sys
lib_path = os.path.abspath(os.path.join('..', 'modules'))
sys.path.append(lib_path)
import fitlorentzian as fl
import getdata as gd

matplotlib.rcParams['font.size'] = 8
matplotlib.rcParams['xtick.major.pad'] = 2
matplotlib.rcParams['ytick.major.pad'] = 2
matplotlib.rcParams['axes.labelpad'] = 2
matplotlib.rcParams['text.usetex'] = True

xlims = (-100, 100)
fisi_susc = (1.7, 1.2)
fisi_bg = (1.7, 1.2)

xlab = 'Probe detuning $\Delta_\mathrm{p} / 2\pi$ (kHz)'
sweep = 'gainsweep'

def col(i):
    N = len(dirs)*1.0
    cm = plt.cm.jet
    return cm((i+1)/N)
red = (0.7, 0.1, 0.1)
dblue = (0.2, 0.2, 0.5)
gray = (0.4, 0.4, 0.4)

opts = {'marker':'.',
        'markersize':2,
        'linestyle':'',
        'rasterized':False,
        }

bpowers, widths, gains, colors, relfreqs, nampses = pickle.load(open('ampl_gainsweep3_shortdat', 'rb'))

fig_namp, ax_namp = plt.subplots(figsize=fisi_susc)
fig_namp.subplots_adjust(0,0,1,1)
ax_namp.set_xlabel(xlab)
ax_namp.set_ylabel('$|S_{11}|^2$ (dB)')
for namps, color in zip(nampses, colors):
        ax_namp.plot(relfreqs, namps, color=color, **opts)
ax_namp.set_xlim(xlims)
ax_namp.set_ylim(0, 45)
ax_namp.locator_params(nbins=5)
figopts = {
        'dpi':400,
        'bbox_inches':'tight' ,
        }
fig_namp.savefig('ampl_{}_vnatraces.pdf'.format(sweep), **figopts)
plt.close(fig_namp)

gbpowers, gwidths, ggains =  pickle.load(open('ampl_gainsweep3_longdat', 'rb'))

inds_not_nan = np.logical_not(np.isnan(ggains))
gbpowers = gbpowers[inds_not_nan]
gwidths = gwidths[inds_not_nan]
ggains = ggains[inds_not_nan]

# fitting bw
gbpowerslin = 10**(gbpowers/10.0)
bpowerslin = 10**(bpowers/10.0)
inds = np.logical_and(gbpowerslin > 0.0, gwidths>0)
a, b = np.polyfit(gbpowerslin[inds], gwidths[inds], 1)
Pthreshold = -b/a
coop = bpowerslin/Pthreshold
gcoop = gbpowerslin/Pthreshold
fcoop = np.linspace(0.45, 1.0, 1001)

fig, ax = plt.subplots(figsize=fisi_bg)
fig.subplots_adjust(0,0,1,1)
ax2 = ax.twinx()
ax, ax2 = ax2, ax # gain on right and bw on left
ax2.set_xlabel(r'Cooperativity $\mathcal{C}$')
ax.set_xlabel(r'Cooperativity $\mathcal{C}$')
ax.set_ylabel(r'Gain $\mathcal{G}$ (dB)  ($\triangle$)')
ax2.set_ylabel(r'Bandwidth (kHz) ($\bigcirc$)')
ms = 12 # marker size
# fitting gaincurve to extract eta
gaincurve = lambda coop, eta: 20*np.log10(np.abs((2*eta-1+coop)/(1-coop)))
fitinds = np.logical_and(gcoop>0.8, ggains<40)
popt, pcov = scipopt.curve_fit(gaincurve, 
        gcoop[fitinds], ggains[fitinds])
feta = popt[0]
ax.plot(fcoop, gaincurve(fcoop, feta), color=dblue)
ax.scatter(coop, gains, color=colors, 
        marker='^', s=ms, zorder=6)
ax.scatter(gcoop, ggains, color=gray, 
        marker='^', s=ms, zorder=5)
ax2.plot(fcoop, b*(1-fcoop), color=dblue)
ax2.scatter(coop, widths, color=colors, 
        marker='o', s=ms, zorder=6)
ax2.scatter(gcoop, gwidths, color=gray, 
        marker='o', s=ms, zorder=5)
ax.locator_params(nbins=5)
ax2.locator_params(nbins=5)
ax.set_xlim(0.49, 1.01)
ax.set_ylim(0, 45) # gain
ax2.set_ylim(0, 70) # bandwidth
fig.savefig('ampl_{}_gainBW.pdf'.format(sweep), **figopts)

