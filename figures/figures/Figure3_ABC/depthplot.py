import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import scipy.optimize as scipopt

matplotlib.rcParams['font.size'] = 8
matplotlib.rcParams['xtick.major.pad'] = 2
matplotlib.rcParams['ytick.major.pad'] = 2
matplotlib.rcParams['axes.labelpad'] = 2
matplotlib.rcParams['text.usetex'] = True

xlabel = r'$-\kappa_\mathrm{om} \, / \, \kappa $'
ylabel = r'$|S_{11}|^2 (\omega=\omega_c)$ (dB)'

# Colors
gray = (0.5, 0.5, 0.5)
lightgray = (0.9, 0.9, 0.9)
palepurple = (0.93, 0.89, 1.0)
paleorange = (1.0, 0.95, 0.9)
dblue = '#0504aa'
dred = '#8f1402'
fitcol = '#014d4e'

def depth(fkap0, fkapdba, fkapc):
    return ((fkap0 + fkapdba - fkapc) / (fkap0 + fkapdba + fkapc))**2

data = np.loadtxt('depth.dat')
ppowers = data[:,0]
fkappas = data[:,1]
depths = data[:,2]
fkappa0 = 120 # fitted resonance width
expfkdba = fkappas - fkappa0
indsred = np.nonzero(expfkdba>0)
indsblue = np.nonzero(expfkdba<0)
indsfit = np.nonzero(abs(expfkdba)>25)

plottedPsblue = [9., 8., 6.,  0., -9.]
plottedPsred = [9., 8., 5.5,  0., -9.]
indsplottedred = np.searchsorted(ppowers[indsred], plottedPsred)
indsplottedblue = np.searchsorted(ppowers[indsblue], plottedPsblue)
def colp(P):
    linP = 10**(P/10.0)
    cm = plt.cm.jet
    return cm(1.0 - linP/10.0)

fitfunc = lambda x, fk0, fk: (depth(fk0, x, fk-fk0)) 
popt, pcov = scipopt.curve_fit(fitfunc,
            expfkdba[indsfit], 10**(depths[indsfit]/10.), p0=[68, 120])
fkappa = popt[1]
fkap0 = popt[0]; fkapc = popt[1]-fkap0

xlims = (-120/fkappa, 80/fkappa)
fkapdbas = fkappa* np.linspace(xlims[0], xlims[1], 501)

fig, ax = plt.subplots(figsize=(3.5, 1.2))
fig.subplots_adjust(0,0,1,1)
ax.set_xlabel(xlabel)
ax.set_ylabel(ylabel)
ax.plot(-fkapdbas /fkappa, 10*np.log10(fitfunc(fkapdbas, *popt)),
        color=fitcol)
ax.scatter(-expfkdba[indsred] /fkappa, depths[indsred], marker='s',
        zorder=5, color=dred)
ax.scatter(-expfkdba[indsred][indsplottedred] /fkappa, 
        depths[indsred][indsplottedred], marker='s',
        zorder=4, color=colp(ppowers[indsred][indsplottedred]))
ax.scatter(-expfkdba[indsblue] /fkappa, depths[indsblue], 
        zorder=5, color=dblue)
ax.scatter(-expfkdba[indsblue][indsplottedblue] /fkappa, 
        depths[indsblue][indsplottedblue], marker='o',
        zorder=4, color=colp(ppowers[indsblue][indsplottedblue]))
ax.locator_params(axis='y', nbins=5)
ax.locator_params(axis='x', nbins=10)
ymin, ymax = (-30, 10)
ax.set_ylim(ymin, ymax)
ax.set_xlim(-xlims[1], -xlims[0])

# lines and annotations
fcrit = (fkap0-fkapc) /fkappa
ax.plot((fcrit, fcrit), (ymin, ymax), '--', color=gray)
xfb = np.linspace(-fkap0 /fkappa, xlims[0], 101)
ax.fill_between(-xfb, -30, 10,
        color=paleorange)
anopts = {
        "xytext":(0,0),
        'textcoords':'offset points',
        'ha':'center',
        'va':'center',
        }
h=3
ax.annotate(('undercoupled\n'+
    r'$\kappa_0+\kappa_\mathrm{om}>\kappa_\mathrm{ex} $'),
        xy=((-xlims[1]+fcrit)/2.0, h), **anopts)
ax.annotate(('overcoupled\n'+
    r'$\kappa_0+\kappa_\mathrm{om}<\kappa_\mathrm{ex} $'),
        xy=((fcrit+fkap0/fkappa)/2.0, h), **anopts)
ax.annotate(('amplification\n'+
    r'$\kappa_0+\kappa_\mathrm{om}<0$'),
        xy=((fkap0/fkappa - xlims[0])/2.0, -10), **anopts)

fig.savefig('depth.pdf', bbox_inches='tight')

