import matplotlib.pyplot as plt
import matplotlib

matplotlib.rcParams['font.size'] = 8
matplotlib.rcParams['xtick.major.pad'] = 2
matplotlib.rcParams['ytick.major.pad'] = 2
matplotlib.rcParams['axes.labelpad'] = 2
matplotlib.rcParams['text.usetex'] = True


def colp(P):
    linP = 10**(P/10.0)
    cm = plt.cm.jet
    return cm(1.0 - linP/10.0)

fig, ax = plt.subplots(figsize=(0.1, 1.))
fig.subplots_adjust(0,0,1,1)

Pmin = 0; Pmax = 10
conv_const = 0.0882186698941 # factor to go from power to coop derived in other script
Cmin = Pmin*conv_const; Cmax = Pmax*conv_const
norm = matplotlib.colors.Normalize(vmin=Cmin, vmax=Cmax)
cb = matplotlib.colorbar.ColorbarBase(ax, cmap='jet_r', norm=norm, 
        orientation='vertical', ticks=[0,0.4,0.8])
ax.set_xlabel(r'$\mathcal{C}$')

fig.savefig('colorbar.pdf', bbox_inches='tight')

