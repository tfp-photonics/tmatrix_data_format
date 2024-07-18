import numpy as np
import os, sys
import h5py
#possibly multiple wavelengths
wavelength  = np.arange(400,850, 50)
Nl = len(wavelength)
lmax = 3
qmax = 2*(lmax*(lmax+1)+lmax)
# dummy 30x30 matrix values repeated for each wavelength
tdata = np.reshape(np.arange(1,901,1) + 1j*np.arange(1,901,1), (qmax,qmax))
tmatrix = np.zeros((Nl,qmax,qmax), dtype="complex")
for i in range(Nl):
  tmatrix[i,:,:] = tdata

print(tmatrix[0,0:3,0:3])

l = np.zeros(qmax)
m = np.zeros(qmax)
s = np.array([b'xxxxxxic']*len(l))
i=0
for li in range(1,lmax+1):
    for mi in range(-li, li+1):
        for si in [b'magnetic', b'electric']:
            l[i] = li
            m[i] = mi
            s[i] = si
            i = i+1

# Saving to HDF5
f = 'ap.tmat.h5';
with h5py.File(f, "w") as f:
    f["vacuum_wavelength"] = np.asarray(wavelength)
    f["vacuum_wavelength"].attrs["unit"] = "nm"
    f["tmatrix"] = tmatrix
    f['modes/l'] = l
    f['modes/m'] = m
    f['modes/polarization'] = s
    emb = f.create_group(f"embedding")
    emb["relative_permittivity"] = 1.33**2
    emb["relative_permeability"] = 1.0
    emb.attrs["name"] = "H2O, Water"
    emb.attrs["keywords"] = "non-dispersive"
    sca_mat = f.create_group("scatterer/material")
    sca_mat.attrs["name"] = 'Au, Gold'
    sca_mat.attrs["keywords"] = "dispersive, plasmonic"
    sca_mat.attrs["reference"] = "Au from Raschke et al 10.1103/PhysRevB.86.235147"
    sca_mat["relative_permittivity"] = np.ones(len(wavelength), dtype = complex)*(-11.4+1.181j)
    sca_mat["relative_permeability"] = 1.0
    sca_gr = f.create_group("scatterer/geometry")
    sca_gr["radiusxy"] = 20.0
    sca_gr["radiusz"] = 40.0
    sca_gr.attrs['unit'] = 'nm'
    sca_gr.attrs['shape'] = 'spheroid'
    sca_gr.attrs['name'] = 'homogeneous spheroid with symmetry axis z'
    mpar = f.create_group("computation/method_parameters")
    #f["computation/analytical_zeros"] = analytical_zeros # not needed in this example
    mpar["Lmax"] = lmax
    mpar["Ntheta"] = 100
    f.create_group("computation/files")
    with open(__file__, "r") as scriptfile:
        f[f"computation/files/{os.path.basename(__file__)}"] = scriptfile.read()
    f["computation"].attrs["method"] = "EBCM, Extended Boundary Condition Method" 
    f["computation"].attrs["description"] = "Computation using SMARTIES, a numerically robust EBCM implementation for spheroids" 
    f["computation"].attrs["software"] = f"SMARTIES=1.1, python={sys.version.split()[0]}, h5py={h5py.__version__}"
    f["computation"].attrs["name"] = "SMARTIES"
    f.attrs['name'] = 'Au prolate spheroid in water'
    f.attrs['keywords'] = 'gold, spheroid, ebcm, passive, reciprocal, czinfinity, mirrorxyz'
    f.attrs['description'] = 'Computation using SMARTIES, a numerically robust EBCM implementation for spheroids'
    f.attrs['storage_format_version'] = 'v1'    

