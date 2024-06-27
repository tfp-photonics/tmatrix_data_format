# demomulti02.py - Scattering spectra for single nanoparticle and planewave excitation.

import numpy as np

import treams
import treams.io
import h5py

# import previously stored T-matrices
finp = "tmatrix_cylinder.h5.h5"
tmat = treams.io.load_hdf5(finp,"nm")

# array for extinction and scattering cross section
n = np.size(tmat)
sig = np.zeros((n,2)) 

# loop over wavenumbers
for i in range(n):
    tm = tmat[i]
    inc = treams.plane_wave([tm.k0, 0, 0], 0, k0=tm.k0, material=tm.material, poltype="parity")
    sca = tm @ inc.expand(tm.basis)
    
    sig[i,0], sig[i,1] = tm.xs(inc)

# save output to file
np.savetxt("spectra02.txt", sig)