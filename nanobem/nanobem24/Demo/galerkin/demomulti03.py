# demomulti02.py - Scattering spectra for coupled nanoparticles.

import numpy as np

import treams
import treams.io
import h5py

# import previously stored T-matrices
finp = "tmatrix_cylinder.tmat.h5"
tmat = treams.io.load_hdf5(finp,"nm")
# sphere positions
shift = 400 * np.array([[ -1, 0, 0 ], [ 1, 0, 0 ]]) 

# array for extinction and scattering cross section
n = np.size(tmat)
sig = np.zeros((n,2)) 

# loop over wavenumbers
for i in range(n):
    tm = tmat[i]
    tm = treams.TMatrix.cluster([ tm, tm ], shift).interaction.solve()
    inc = treams.plane_wave([0, 0, tm.k0], [ 1, 0, 0 ], k0=tm.k0, 
                                  material=tm.material, poltype="parity")
    sca = tm @ inc.expand(tm.basis)
    
    sig[i,0], sig[i,1] = tm.xs(inc)


# save output
np.savetxt("spectra03.txt", sig)