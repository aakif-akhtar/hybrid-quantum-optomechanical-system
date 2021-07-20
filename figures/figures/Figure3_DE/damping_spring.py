import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import os.path as path
import scipy.optimize as scipopt

matplotlib.rcParams['font.size'] = 8
matplotlib.rcParams['xtick.major.pad'] = 2
matplotlib.rcParams['ytick.major.pad'] = 2
matplotlib.rcParams['axes.labelpad'] = 2
matplotlib.rcParams['text.usetex'] = True

xlab = 'Pump detuning $\Delta \, / \, 2 \pi$ (MHz)'
xl1 = (-5.95, -4.7)
xl2 = (4.7, 5.95)
pxl1 = (-5.95, -4.7)
pxl2 = (4.7, 5.95)

def getsourceP(dirname, srcst):
    basename = path.basename(dirname)
    with open(dirname + '/' + basename + '.set', 'r') as setfile:
        for line in setfile:
            if srcst in line:
                next(setfile); next(setfile)
                power = float(next(setfile).split()[1])
    return power

def brokenaxes(fig, ax, ax2):
    # reduce ticks
    ax.locator_params(axis = 'x', tight = True, nbins = 5)
    ax2.locator_params(axis = 'x', tight = True, nbins = 5)
    # The following arbitrarily moves xlabel y and x
    # Optimally, it should only move according to x
    ax.xaxis.set_label_coords(0.5, -0.135, transform=fig.transFigure)
    # hide the spines between ax and ax2
    ax.spines['right'].set_visible(False)
    ax2.spines['left'].set_visible(False)
    ax.yaxis.tick_left()
    ax.tick_params(labelright='off')
    ax2.yaxis.tick_right()
    ax2.tick_params(labelright='off')
    d = .015 # how big to make the diagonal lines in axes coordinates
    kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
    # diagonal lines
    ax.plot((1-d,1+d), (-d,+d), **kwargs)
    ax.plot((1-d,1+d),(1-d,1+d), **kwargs)
    kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
    ax2.plot((-d,+d), (1-d,1+d), **kwargs)
    ax2.plot((-d,+d), (-d,+d), **kwargs)
    fig.subplots_adjust(wspace=0.05)

def fkappa(fDelta, fkappa_0, A, fGammam, fOmegam):
    deltafkappa = A * (
            1.0/ (1.0 + ((fDelta+fOmegam)/(fGammam/2.0))**2)
            -
            1.0/ (1.0 + ((fDelta-fOmegam)/(fGammam/2.0))**2)
            )
    return fkappa_0 + deltafkappa

def fres(fDelta, fres_0, B, fGammam, fOmegam):
    deltafres = - B * (
            (fDelta + fOmegam) / (1 + ((fDelta+fOmegam)/(fGammam/2.0))**2)
            -
            (fDelta - fOmegam) / (1 + ((fDelta-fOmegam)/(fGammam/2.0))**2)
            )
    return fres_0 + deltafres

dirs = ['031422_mechanical_damping_deltasweep_probe_5dBm']

for d in dirs:
    basename = path.basename(d)
    data = np.loadtxt(d + '/' + basename + '.dat')

    powr = getsourceP(d, 'src2')
    plotsuf = basename.split('_')[0] + '_{}dBm'.format(powr)

    freqs = data[:,0]
    relfreqs = data[:,1]
    width = data[:,2]
    fshift = data[:,3]
    fshift0 = fshift[0]
    fshift = fshift-fshift0
    fshift = 1e-3*fshift # kHz
    R2 = data[:,4]

    xlim = [-8, 8]
    inds = np.nonzero(np.logical_and(relfreqs>xlim[0], relfreqs<xlim[1]))
    freqs = freqs[inds]
    relfreqs = relfreqs[inds]
    width = width[inds]
    fshift = fshift[inds]
    R2 = R2[inds]

    red = (0.7, 0.1, 0.1)
    dblue = (0.2, 0.2, 0.5)
    opts = {
            'marker':'o',
            's':10,
            'zorder':1,
            'color': red,
            }
    fitopts = {'marker':'',
            'linestyle':'-',
            'color': dblue,
            }
    fisi = (2.5, 1.2)

    figopts = {
            'dpi':400,
            'bbox_inches':'tight',
            }

    # fitting
    fitfreqs1 = np.linspace(-5.82, -4.81, 1001)
    fitfreqs2 = np.linspace(4.78, 5.78, 1001)
    finds1 = np.nonzero(np.logical_and(relfreqs>xl1[0], relfreqs<xl1[1]))
    finds2 = np.nonzero(np.logical_and(relfreqs>xl2[0], relfreqs<xl2[1]))
    # mechanical damping
    fkappa00 = 115; fOmegam0 = 5.33; fGammam0 = 0.5; A0 = 30;
    p0 = [fkappa00, A0, fGammam0, fOmegam0]
    pdamp1, covdamp1 = scipopt.curve_fit(fkappa, relfreqs[finds1], width[finds1],
            p0=p0)
    pdamp2, covdamp2 = scipopt.curve_fit(fkappa, relfreqs[finds2], width[finds2],
            p0=p0)
    # mechanical spring
    B0 =  0.9* pdamp2[1] / pdamp2[2] # B = A / Gammam
    fres0 = 10;
    p0_s = [fres0, B0, pdamp2[2], pdamp2[3]]
    try:
        pspring1, covspring1 = scipopt.curve_fit(fres, relfreqs[finds1], fshift[finds1],
                p0=p0_s)
        pspring2, covspring2 = scipopt.curve_fit(fres, relfreqs[finds2], fshift[finds2],
                p0=p0_s)
    except Exception as e:
        print 'exception'
        print e
        pspring1 = p0_s
        pspring2 = p0_s


    # Broken axes
    fig, (ax, ax2) = plt.subplots(1,2, sharey=True, figsize=fisi)
    fig.subplots_adjust(0,0,1,1)
    ax.set_xlim(pxl1)
    ax2.set_xlim(pxl2)
    brokenaxes(fig, ax, ax2)
    ax.set_xlabel(xlab)
    ax.set_ylabel(r'$(\kappa + \kappa_\mathrm{om}) \, / \, 2\pi$ (kHz)')
    ax.scatter(relfreqs, width, **opts)
    ax.plot(fitfreqs1, fkappa(fitfreqs1, *pdamp1), **fitopts)
    ax2.scatter(relfreqs, width, **opts)
    ax2.plot(fitfreqs2, fkappa(fitfreqs2, *pdamp2), **fitopts)
    ax.locator_params(axis = 'y', tight = True, nbins = 5)
    ax.locator_params(axis = 'x', tight = True, nbins = 3)
    ax2.locator_params(axis = 'x', tight = True, nbins = 3)
    adjbot = 0.25
    fig.savefig(plotsuf + 'brokenaxes_width.pdf', **figopts)
    plt.close(fig)

    # Broken axes
    fig, (ax, ax2) = plt.subplots(1,2, sharey=True, figsize=fisi)
    fig.subplots_adjust(0,0,1,1)
    ax.set_xlim(pxl1)
    ax2.set_xlim(pxl2)
    brokenaxes(fig, ax, ax2)
    ax.set_xlabel(xlab)
    ax.set_ylabel(r'$\delta \omega_\mathrm{om} \, / \, 2\pi$ (kHz)')
    ax.scatter(relfreqs, fshift, **opts)
    ax.plot(fitfreqs1, fres(fitfreqs1, *pspring1), **fitopts)
    ax2.scatter(relfreqs, fshift, **opts)
    ax2.plot(fitfreqs2, fres(fitfreqs2, *pspring2), **fitopts)
    ax.set_ylim(-10, 20)
    ax.locator_params(axis = 'y', tight = True, nbins = 4)
    ax.locator_params(axis = 'x', tight = True, nbins = 3)
    ax2.locator_params(axis = 'x', tight = True, nbins = 3)
    fig.savefig(plotsuf + 'brokenaxes_fshift.pdf', **figopts)
    plt.close(fig)

