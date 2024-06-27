# demomulti01.py - Averaged scattering spectra for single nanoparticle.

import numpy as np

import treams
import treams.io
import h5py

# import previously stored T-matrices
finp = "tmatrix_cylinder.h5.h5"
tmat = treams.io.load_hdf5(finp,"nm")

# array for average extinction and scattering cross section
n = np.size(tmat)
sig = np.zeros((n,2)) 

# loop over wavenumbers
for i in range(n):
    sig[i,0] = tmat[i].xs_sca_avg
    sig[i,1] = tmat[i].xs_ext_avg

# save output to file
np.savetxt("spectra01.txt", sig)