import matplotlib
import matplotlib.pyplot as plt
import os
from math import sqrt

def smallsettings():
    matplotlib.rcParams['font.size'] = 12
    matplotlib.rcParams['xtick.major.pad'] = 2
    matplotlib.rcParams['ytick.major.pad'] = 2
    matplotlib.rcParams['axes.labelpad'] = 2
    matplotlib.rcParams['text.usetex'] = True
    matplotlib.rcParams["pdf.fonttype"] = 42

figopts = {
        'figsize':(1.7, 1.2),
        }

def ensuredir(fname):
    di = os.path.dirname(fname)
    if not di == '':
        if not os.path.exists(di):
            os.makedirs(di)

def graph(xx,yy, name, xlab='x', ylab='y', legends=None,
        logx=False, logy=False, colors=None,
        fit=None, yrange=None, xrange=None, 
        points=False, ps=12, hline=None, vlines=None, small=False,
        fbetween=None, otherfunc=None, rast=True, dpi=400,
        errors=None, saveopts={}, fisi=(1.7, 1.2)):
    """
    Small wrapper for simple plots with matplotlib
    """
    if small:
        smallsettings()
    fig, ax = plt.subplots(figsize=fisi)
    fig.subplots_adjust(0,0,1,1)
    plotopts = {}
    if rast:
        plotopts['rasterized']=True
    if points:
        plotopts['linestyle'] = 'None'
        plotopts['marker'] = 'o'
        plotopts['markersize'] = sqrt(ps)
        plotopts['markeredgecolor'] = 'none'
    if hline is not None:
        lopt = 'k--'
        if xrange is not None:
            ax.plot((xrange[0], xrange[1]), (hline, hline), lopt)
        else:
            ax.plot((min(xx[0]), max(xx[0])), (hline, hline), lopt)
    if vlines is not None:
        for v in vlines:
            lopt = 'k--'
            if yrange is not None:
                ax.plot((v, v), (yrange[0], yrange[1]), lopt)
            else:
                ax.plot((v, v), (min(yy[0]), max(yy[0])), lopt)
    ps = []
    if colors is not None:
        ax.set_color_cycle(colors)
    if fit is not None:
        ax.plot(fit[0], fit[1])
    if errors is None:
        for x, y in zip(xx, yy):
            ps.append(ax.plot(x, y, **plotopts)[0])
    else:
        for x, y, e in zip(xx, yy, errors):
                ps.append(ax.errorbar(x, y, yerr=e, **plotopts))
    if legends is not None:
        plt.legend(ps, legends, frameon=False)
    if fbetween is not None:
        ax.fill_between(fbetween['x'], fbetween['y0'], fbetween['y1'], 
                facecolor=fbetween['col'], 
                alpha=0.5, interpolate=True, linewidth=0.0,)
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    if small:
        ax.locator_params(nbins=5)
    if logx:
        ax.set_xscale('log', nonposx='clip')
    if logy:
        ax.set_yscale('log', nonposx='clip')
    if yrange is not None:
        ax.set_ylim(yrange)
    if xrange is not None:
        ax.set_xlim(xrange)
    if otherfunc is not None:
        otherfunc(fig, ax)
    ensuredir(name)
    if small:
        fig.savefig(name, bbox_inches='tight', dpi=dpi, **saveopts)
    else:
        fig.savefig(name, **saveopts)
    plt.close(fig)

def colormap(Z, gname, ran=None,
        xlab='x', ylab='y', clab='c',
        cmap='RdBu_r', lims=None, vslices=[], hslices=[], fisi=(3,2),
        hide_cbar=False, hide_ylabs=False, vlims=None, otherfunc=None,
        small=False):
    if small:
        smallsettings()
    options = {
            'aspect':'auto',
            'origin':'lower',
            'interpolation':'nearest', # nearest is what to use for pdf or svg
            }
    if vlims is not None:
        options['vmin'] = vlims[0]
        options['vmax'] = vlims[1]
    xmin, xmax, ymin, ymax = None, None, None, None
    if ran is not None:
        options['extent'] = ran
    fig, ax = plt.subplots(figsize=fisi)
    fig.subplots_adjust(0,0,1,1)
    plotamp = ax.imshow(Z, **options) 
    plotamp.set_cmap(cmap)
    cbar = plt.colorbar(plotamp, shrink=1)
    cbar.set_label(clab)
    if hide_cbar:
        cbar.remove()
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    if hide_ylabs:
        ax.yaxis.set_ticklabels([])
        ax.set_ylabel('')
    if lims:
        ax.set_xlim((lims[0], lims[1]))
        ax.set_ylim((lims[2], lims[3]))
    ax.locator_params(nbins=5)
    for vslice in vslices:
        ax.plot((vslice, vslice), ax.get_ylim(),
                '--', color='grey', lw=0.5,
                )
    for hslice in hslices:
        ax.plot(ax.get_xlim(), (hslice, hslice),
                '--', color='grey', lw=0.5,
                )
    if otherfunc is not None:
        otherfunc(fig, ax)
    ensuredir(gname)
    fig.savefig(gname, bbox_inches='tight', dpi=600)
    plt.close(fig)

