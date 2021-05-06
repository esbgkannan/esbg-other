import pickle
import numpy as np
from pyemma.thermo import estimate_umbrella_sampling as estimate_us

K = 178
kB = 1.381e-23 * 6.022e23 / 1000.0 # Boltzmann constant in kJ/mol/K
temperature = 310 
kt = kB*temperature
beta = 1/kt

omit = 500 # omit the first 1000 dataset
us_trajs = []
with open('us_data1_cv.txt') as cvfiles:
#with open('us_data2_cv.txt') as cvfiles:
#with open('us_data3_cv.txt') as cvfiles:
    for line in cvfiles.readlines():
        this_cv = np.loadtxt(line.strip())
        us_trajs.append(np.array(this_cv[omit:,1:3], order='C'))

us_centers = []
us_force_constants = []
with open('spring_order.txt') as params:
    for line in params.readlines():
        c1, c2, k1, k2 = line.split()
        us_centers.append(np.array([float(c1), float(c2)]))
        us_force_constants.append(np.array([float(k1), float(k2)]))

cv = np.concatenate(us_trajs)
pd_kn = np.reshape([cv[:,0]], (K, -1))
ke_kn = np.reshape([cv[:,1]], (K, -1))

pd_min = 0.15
pd_max = 3.0
ke_min = -1.9
ke_max = 2.1
nbins_per_cv = 50
dx1 = (pd_max - pd_min) / float(nbins_per_cv)
dx2 = (ke_max - ke_min) / float(nbins_per_cv)
bin_kn = np.ones(pd_kn.shape, np.int16) * (-1) # initialize with -1 so that we can easily identify points that are not assigned to a bin
nbins = 0
bin_counts = list()
bin_centers = list()
not_selected_bin_centers = list()
for i in range(nbins_per_cv):
    for j in range(nbins_per_cv):
        pd = pd_min + dx1 * (i + 0.5)
        ke = ke_min + dx2 * (j + 0.5)
        in_bin = (pd-dx1/2 <= pd_kn) & (pd_kn <= pd+dx1/2) & (ke-dx2/2 <= ke_kn) & (ke_kn <= ke+dx2/2)
        bin_count = in_bin.sum()
        if (bin_count > 0): 
            bin_centers.append((pd, ke))
            bin_counts.append(bin_count)
            bin_kn[in_bin] = nbins
            nbins += 1
        else:
            not_selected_bin_centers.append((pd, ke))

us_dtrajs = []
for i in range(bin_kn.shape[0]):
    us_dtrajs.append(bin_kn[i,:])

wham = estimate_us(us_trajs, us_dtrajs, us_centers, us_force_constants, maxiter=100000, maxerr=1.0E-15, save_convergence_info=1, estimator='wham', kT=kt)

f = open('wham.p', 'w')
pickle.dump(wham, f)
f.close()
output = ''
for i in range(nbins):
    output+='%8d %6.3f %6.3f %8d %10.3f\n' % (i, bin_centers[i][0], bin_centers[i][1], bin_counts[i], wham.free_energies[i])
for ct in not_selected_bin_centers:
    i += 1
    output+='%8d %6.3f %6.3f %8d inf\n' % (i, ct[0], ct[1], 0)
with open('wham_pmf.txt', 'w') as out:
    out.write(output)
