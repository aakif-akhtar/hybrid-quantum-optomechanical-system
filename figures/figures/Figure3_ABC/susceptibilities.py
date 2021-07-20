import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import os, sys
lib_path = os.path.abspath(os.path.join('..', 'modules'))
sys.path.append(lib_path)
import fitlorentzian as fl

matplotlib.rcParams['font.size'] = 8
matplotlib.rcParams['xtick.major.pad'] = 2
matplotlib.rcParams['ytick.major.pad'] = 2
matplotlib.rcParams['axes.labelpad'] = 2
matplotlib.rcParams['text.usetex'] = True

xlab = 'Probe detuning $\Delta_\mathrm{p} / 2\pi$ (kHz)'

swpdirs = ['134757_vna_Pump power (dBm)', 
        '143732_vna_Pump power (dBm)']

for swpd in swpdirs:
    sweep = swpd
    fulldata = np.loadtxt(swpd + '/' + swpd + '.dat')
    Npowers = len(np.unique(fulldata[:,0]))

    if '13' in swpd:
        amp=True
        plottedPs = [9., 8., 6.,  0., -9.]
    else:
        amp=False
        plottedPs = [9., 8., 5.5,  0., -9.]

    def col(i):
        N = Npowers *1.0
        cm = plt.cm.jet
        return cm(1.0-i/N)
    def colp(P):
        linP = 10**(P/10.0)
        cm = plt.cm.jet
        return cm(1.0 - linP/10.0)
    red = (0.7, 0.1, 0.1)
    dblue = (0.2, 0.2, 0.5)
    opts = {
            'marker':'.',
            'markersize':2,
            'linestyle':'',
            'rasterized':False,
            }

    figopts = {
            'figsize':(1.5, 1.2),
            }
    ylab = '$|S_{11}|^2$ (dB)'

    fig_namp, ax_namp = plt.subplots(**figopts)
    fig_namp.subplots_adjust(0,0,1,1)
    ax_namp.set_xlabel(xlab)
    if not amp:
        ax_namp.set_ylabel(ylab)

    centfreq = None
    for i, pdat in enumerate(np.split(fulldata, Npowers)):
        if len(np.unique(pdat[:,0]))>1: print 'Problem'
        bpump = pdat[0,0]
        if bpump <= 10 and bpump > -10:

            freqs = pdat[:,1]
            logamps = pdat[:,2]
            phase = pdat[:,4]

            if centfreq is None:
                centfreq0 = freqs[len(freqs)//2]
                relfreqs0 = 1e6* (freqs - centfreq0) 
                normpar = (fl.fitting(relfreqs0, 10**(logamps/10.0)))
                centfreq = centfreq0 +  normpar[2]*1e-6

            relfreqs = 1e6* (freqs - centfreq) 

            bg = 10*np.log10(normpar[0])
            namps = logamps - bg

            if bpump in plottedPs:
                ax_namp.plot(relfreqs, namps, color=colp(bpump), **opts)

    if amp:
        ax_namp.yaxis.set_ticklabels([])
    ax_namp.set_xlim((-150, 150))
    ax_namp.set_ylim((-25, 5))
    ax_namp.locator_params(nbins=4)
    sfigopts = {
            'dpi':400,
            'bbox_inches':'tight',
            }
    fig_namp.savefig(sweep + '_amplitude_normalized.pdf', **sfigopts)


