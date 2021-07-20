import os, sys
lib_path = os.path.abspath(os.path.join('..', 'modules'))
sys.path.append(lib_path)
from graph import graph
import pickle
import numpy as np

def gain_cf(res):
    bg = np.abs(res['Renormalisation'])
    amplitude = (2*res['r0'] - bg)**2 - bg**2
    return 10*np.log10(amplitude / bg**2)

pumpsP, linpumpsP, gainfuncs, fkappa_effs, coops, cooppars, coops_err, fres = pickle.load(
	open('Sweep_gaindata', 'rb'))

# Gain curve
gains = np.array( map(gain_cf, gainfuncs) )
inds = coops>0.5

neffdatname = 'Sweep_neffdata'
neffs, N_omcs = pickle.load(open(neffdatname, 'rb'))

# fitting added noise
def addednoise(C, eta, neff, n0):
    dem = 1./ (1-C-2*eta)**2
    term1 = 4*C*eta*(neff + 0.5)
    term2 = 4 * eta * (1-eta) * (n0 + 0.5)
    return (term1 + term2) * dem

# figure noise graph
xrange=[20, 45]
feta = 150./197. # eta value from fitting resonances
def sql_line(fig, ax):
    hline=addednoise(1, feta, 0,0)
    ax.plot((xrange[0], xrange[1]), (hline, hline), '--', color='#ec2d01')
    ax.text(43, hline+0.05, r'\textit{device quantum limit}',
            va='bottom', ha='right', fontstyle='italic',
            )

graph([gains[inds]], [N_omcs[inds]],
    'Noise.pdf',
    xlab='Gain (dB)', ylab='$\mathcal{N}(\omega_c)$ (photons/(s Hz))', 
    xrange=xrange, yrange=[0,3],
    colors=['#8c000f'],
    otherfunc=sql_line,
    points=True, small=True, rast=False) 

