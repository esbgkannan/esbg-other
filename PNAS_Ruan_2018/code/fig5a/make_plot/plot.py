#!/usr/bin/python2.7

"""python script to generate fig5a

Make sure that files in chelix_rmsd and ke_distance 
are generated prior to executing the script.
"""

import numpy as _np
import matplotlib.pylab as plt
from matplotlib import cm

def get_x_y_z(xall, yall, nbins=1000, weights=None, kT=1.0, range=[[0,20], [0,3]]):
    z, xedge, yedge = _np.histogram2d(xall, yall, bins=nbins, weights=weights, range=range)
    x = 0.5*(xedge[:-1] + xedge[1:])
    y = 0.5*(yedge[:-1] + yedge[1:])
    # avoid zeros
    zmin_nonzero = _np.min(z[_np.where(z > 0)])
    z = _np.maximum(z, zmin_nonzero)/_np.sum(z)
    # compute free energies
    F = -kT * _np.log(z)
    return x, y, F.T-F.T.min()
    #return x, y, -z.T

def get_data(ke_file, chelix_file):
    ke = [float(i) for i in open(ke_file, 'r').readlines()]
    chelix = [float(i) for i in open(chelix_file, 'r').readlines()]
    x, y, z = get_x_y_z(ke, chelix)
    return x, y, z

# obtain mesh data
wt_x, wt_y, wt_z = get_data('../ke_distance/WT_EGFR_ke.txt',
                            '../chelix_rmsd/WT_EGFR_chelix_rmsd.txt')

g_x, g_y, g_z = get_data('../ke_distance/D770_N771insG_ke.txt',
                         '../chelix_rmsd/D770_N771insG_chelix_rmsd.txt')

gy_x, gy_y, gy_z = get_data('../ke_distance/D770_GY_ke.txt',
                            '../chelix_rmsd/D770_GY_chelix_rmsd.txt')

n_x, n_y, n_z = get_data('../ke_distance/N771_P772insN_ke.txt',
                         '../chelix_rmsd/N771_P772insN_chelix_rmsd.txt')


# set up figure
fig = plt.figure()
c_map = cm.spectral
min_v = min([wt_z.min(), n_z.min(), gy_z.min(), g_z.min()])
max_v = min([wt_z.max(), n_z.max(), gy_z.max(), g_z.max()])
# set the maximum value to 10 for consistensy
wt_z[_np.where(wt_z == wt_z.max())] = 10
n_z[_np.where(n_z == n_z.max())] = 10
gy_z[_np.where(gy_z == gy_z.max())] = 10
g_z[_np.where(g_z == g_z.max())] = 10
norm = cm.colors.Normalize(vmax=max_v, vmin=min_v)

plt.subplot(2, 2, 1)
wt_ax = plt.gca()
wt_ax.set_xlim(left=0.0, right=21)
wt_ax.set_ylim(top=3, bottom=0)
wt_cs = plt.contourf(wt_x*10, wt_y*10, wt_z, 100, cmap=c_map, norm=norm, zorder=10)
wt_ax.spines['top'].set_visible(True)
wt_ax.spines['right'].set_visible(True)
wt_ax.spines['left'].set_visible(True)
wt_ax.spines['bottom'].set_visible(True)
wt_ax.set_title('WT EGFR')

plt.subplot(2, 2, 4)
n_ax = plt.gca()
n_ax.set_xlim(left=0.0, right=21)
n_ax.set_ylim(top=3, bottom=0)
n_cs = plt.contourf(n_x*10, n_y*10, n_z, 100, cmap=c_map, norm=norm, zorder=10)
n_ax.spines['top'].set_visible(True)
n_ax.spines['right'].set_visible(True)
n_ax.spines['left'].set_visible(True)
n_ax.spines['bottom'].set_visible(True)
n_ax.set_title('N771_P772insN')
#cbar = plt.colorbar(n_cs)

plt.subplot(2, 2, 3)
gy_ax = plt.gca()
gy_ax.set_xlim(left=0.0, right=21)
gy_ax.set_ylim(top=3, bottom=0)
gy_cs = plt.contourf(gy_x*10, gy_y*10, gy_z, 100, cmap=c_map, norm=norm, zorder=10)
gy_ax.spines['top'].set_visible(True)
gy_ax.spines['right'].set_visible(True)
gy_ax.spines['left'].set_visible(True)
gy_ax.spines['bottom'].set_visible(True)
gy_ax.set_title('D770>GY')

plt.subplot(2, 2, 2)
g_ax = plt.gca()
g_ax.set_xlim(left=0.0, right=21)
g_ax.set_ylim(top=3, bottom=0)
g_cs = plt.contourf(g_x*10, g_y*10, g_z, 100, cmap=c_map, norm=norm, zorder=10)
g_ax.spines['top'].set_visible(True)
g_ax.spines['right'].set_visible(True)
g_ax.spines['left'].set_visible(True)
g_ax.spines['bottom'].set_visible(True)
g_ax.set_title('D770_N771insG')

#plt.subplot(1,1,1)
#wt_ax = plt.gca()
#wt_ax.set_xlim(left=0.0, right=2.1)
#wt_ax.set_ylim(top=0.3, bottom=0)
#wt_cs = plt.contourf(wt_x, wt_y, wt_z, 100, cmap=c_map, norm=norm, zorder=10)
#plt.colorbar(wt_cs, ticks=[0,1,2,3,4,5,6])
plt.tight_layout()
#plt.savefig('colorbar.pdf', format="pdf")
plt.savefig('fep.pdf', format="pdf")
