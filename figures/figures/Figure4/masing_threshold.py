import numpy as np
import scipy
import matplotlib
import matplotlib.pyplot as plt
import scipy.optimize as scipopt
import pickle
import os, sys
lib_path = os.path.abspath(os.path.join('..', 'modules'))
sys.path.append(lib_path)
import fitlorentzian as fl
import calibpower as cp
import getdata as gd
from graph import graph

matplotlib.rcParams['font.size'] = 8
matplotlib.rcParams['xtick.major.pad'] = 2
matplotlib.rcParams['ytick.major.pad'] = 2
matplotlib.rcParams['axes.labelpad'] = 2
matplotlib.rcParams['text.usetex'] = True

col_amp =  (0.7, 0.2, 0.7)
col_mas =  (0.2, 0.7, 0.7)
blue = '#0504aa'
palepurple = (0.93, 0.89, 1.0)

cmap = 'viridis'
fisi = (2.5, 1.2)
xlength= 1.9
fisi_map = (2.5, 1.2)

xlims = (-200, 200)

BG = -129.142 # HEMT background
calopts = {'HEMT_bg':BG} 

emitlabel = r'$\bar{S}_{aa}(\omega)$ (photons / (s Hz))'
xlab = 'Frequency offset from $\omega_\mathrm{c}/2\pi$ (kHz)'

def col(i):
    N = len(dirs)*1.0
    cm = plt.cm.jet
    return cm(1-i/N)
red = (0.7, 0.1, 0.1)
gray = (0.9, 0.9, 0.9)

fig, ax = plt.subplots(figsize=fisi_map)
fig.subplots_adjust(0,0,1,1)
ax.set_xlabel(xlab)
ax.set_ylabel(r'Pump power $P\, / \, P_\mathrm{th}$') 

i_amp = 24
i_mas = 25

masingdname = 'masingdat'
freqs, powers, widths, peaks, powerdensities, P_amp, fswpowers_amp, P_mas, fswpowers_mas = pickle.load(open(masingdname, 'rb'))

# fitting bw
thpow = 7.75 # threshold power
thpowlin = 10**(thpow/10.0)
bpowerslin = 10**(powers/10.0)
plotinds = bpowerslin < 10**(thpow/10.0)
inds = np.logical_and(
        bpowerslin < thpowlin, 
        bpowerslin > 0.8 * thpowlin)
a, b = np.polyfit(bpowerslin[inds], widths[inds], 1)
Pthreshold = -b/a
coop = bpowerslin/Pthreshold

relfreqs = 1e6*(freqs - freqs[len(freqs)//2])
extents = (relfreqs.min(), relfreqs.max(), min(coop), max(coop))
Z0 = cp.FSWpower_to_Saa(np.array(powerdensities), **calopts)
Z0[Z0<=0] = 1e-3
Z = np.log10(Z0)
im = ax.imshow(Z, aspect='auto', cmap=cmap, extent=extents,
        origin='lower', interpolation='nearest', vmin=0,
        rasterized=False)
cbaxes = fig.add_axes([(2.4)/xlength, 0.0, 0.1/xlength, 1])
cbar = plt.colorbar(im, cax=cbaxes, ticks=[0, 5, 10]) 
cbar.ax.yaxis.tick_left()
cbar.ax.yaxis.set_label_position('left')
cbar.ax.set_yticklabels(['$10^0$', '$10^5$', '$10^{10}$'])
cbar.set_label(emitlabel)
ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
figopts = {
        'dpi':1200,
        'bbox_inches':'tight' ,
        }
c_amp = 10**(P_amp/10)/Pthreshold
c_mas = 10**(P_mas/10)/Pthreshold
ax.plot((relfreqs.min(), relfreqs.max()), (c_amp, c_amp),
        '--', color=col_amp)
ax.plot((relfreqs.min(), relfreqs.max()), (c_mas, c_mas),
        '--', color=col_mas)
ax.set_xlim(xlims)
ax.set_ylim((min(coop), max(coop)))
ax.locator_params(axis='x',nbins=5)
ax.locator_params(axis='y',nbins=8)
fig.savefig('masing_map.pdf', **figopts)
plt.close(fig)

# reference (zero power)
reffreqs, reffswpowers, refparams = gd.getFSWdata(
        '100736_fsw_trace')
col_ref = '#0a437a'

fig, ax = plt.subplots(figsize=fisi)
fig.subplots_adjust(0,0,1,1)
ax.set_xlabel(xlab)
ax.set_ylabel(emitlabel)
opts = {'linewidth':2, 'rasterized':False}
ax.plot(relfreqs, cp.FSWpower_to_Saa(fswpowers_mas, **calopts),
        '-', color=col_mas, **opts)
ax.plot(relfreqs, cp.FSWpower_to_Saa(fswpowers_amp, **calopts),
        '-', color=col_amp, **opts)
# integrate
phots_amp = scipy.integrate.simps(cp.FSWpower_to_Saa(fswpowers_amp, **calopts), relfreqs*1e3)
phots_mas = scipy.integrate.simps(cp.FSWpower_to_Saa(fswpowers_mas, **calopts), relfreqs*1e3)
# fit lorentzian under threshold
bpar = (fl.fitting(1e9*freqs, 
    cp.FSWpower_to_Saa(fswpowers_amp, **calopts)))
rate = bpar[1]; fkap_eff=bpar[3]
neff = rate * (2*np.pi*fkap_eff) / (c_amp * 2*np.pi*150e3 * 2*np.pi*200e3) - 0.5
ax.plot(relfreqs, fl.lorentzian(1e9*freqs, bpar),
        '-', color='#6258c4')

# reference level
ax.plot(relfreqs, cp.FSWpower_to_Saa(reffswpowers, **calopts),
        '-', color=col_ref, **opts)

ax.locator_params(nbins=4)
ax.set_yscale('log')
ax.set_xlim(xlims)
ax.set_ylim(1, 2*10**10)
ax.set_yticks(10**np.array(range(11)))
ax.set_yticklabels( ['$10^0$'] + ['']*4 +
        ['$10^5$'] + ['']*4 + ['$10^{10}$'])
insax = fig.add_axes([0.63, 0.5, 0.3, 0.4])
insax.plot(coop[plotinds], widths[plotinds], 
        '.', markersize=4, color=blue)
endp = 0.3; maxup = 70
insax.plot((1,1),(0, maxup), ':', color='k')
insax.fill_between((1, 1.+endp), (maxup, maxup), color=palepurple)
insax.tick_params('both', length=2, width=0.5, which='major')
small = 6
insax.annotate('masing', xy=(1+endp/2.0, maxup/2.0), size=small,
    ha='center', va='center', rotation='vertical')
insax.set_xlim(0.5, 1+endp)
insax.set_ylim(0, maxup)

matplotlib.rcParams['xtick.major.pad'] = 0.5
matplotlib.rcParams['ytick.major.pad'] = 0.5
matplotlib.rcParams['axes.labelpad'] = 0.5
insax.set_xlabel(r'Cooperativity $\mathcal{C}$')
insax.set_ylabel('Width (kHz)')
insax.locator_params(nbins=3)
for item in ([insax.xaxis.label, insax.yaxis.label] +
        insax.get_xticklabels() + insax.get_yticklabels()):
    item.set_fontsize(small)
fig.savefig('masing_examples.pdf', **figopts)
plt.close(fig)


