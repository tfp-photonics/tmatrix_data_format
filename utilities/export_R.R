## ----eval=FALSE---------------------------------------------------------------
library(rhdf5) # note: dev version to support complex
library(uuid) # uuid
library(glue) # string interpolation
library(purrr) # mapping functions

## dummy data

# possibly multiple wavelengths
wavelength <- seq(400, 800, by=50)
Nl <- length(wavelength)

Lmax <- 3
qmax <- 2*(Lmax*(Lmax+1)+Lmax) # T-matrix size

# dummy 30x30 matrix values for each wavelength
# note the byrow due to HDF5 expecting
# row-major ordering vs R's default column-major
tdata <- matrix(1:qmax^2 + 1i*(1:qmax^2), qmax, qmax, byrow=TRUE)

tmatrix <- array(NA_complex_, c(Nl,qmax,qmax))
for(i in seq_len(Nl))
  tmatrix[i,,] <- tdata

tmatrix[1,1:3,1:3]

modes <- list(l = rep(NA_integer_, qmax),
              m = rep(NA_integer_, qmax),
              polarization = rep(NA_character_, qmax))

i <- 1
for (li in 1:Lmax){
  for (mi in -li:li){
    for (si in c("magnetic", "electric")){
      modes$l[i] <- li
      modes$m[i] <- mi
      modes$polarization[i] <- si
      i <-  i+1
    }
  }
}

embedding <- list('relative_permeability' = 1.0, 
                  'relative_permittivity' = 1.33^2)
scatterer <- list(material = list('relative_permeability' = 1.0, 
                                  'relative_permittivity' =
                                    rep(-11.4+1.181i,Nl)),
                  geometry = list('radiusxy' = 20.0, 
                                  'radiusz' = 40.0))

computation <- list(method_parameters = list('Lmax' = Lmax, 
                                             'Ntheta' = 100),
                    files = list(script = paste(readLines('export_R.R'), collapse = "\n")))



## ----eval=FALSE---------------------------------------------------------------
library(rhdf5) # note: requires dev version to support complex arrays
# cf https://support.bioconductor.org/p/9156305/#9156408
# install.packages("remotes")
# remotes::install_github("grimbough/rhdf5@devel")
# may need 'crypto' library from openSSL to compile
# e.g. via brew install openssl on macos


## ----eval=FALSE---------------------------------------------------------------
f <- 'ar.tmat.h5'
software = sprintf("SMARTIES=1.1, R=%s, rhdf5=%s", paste0(R.version$major, R.version$minor), packageVersion("rhdf5"))

unlink(f) # delete previous file if it exists
h5createFile(f)
h5closeAll() # in case connections open

h5write(wavelength, file=f, name="/vacuum_wavelength")

h5write(tmatrix, file=f, name="/tmatrix", native=TRUE) # store rowwise

h5write(modes, file=f, name='/modes')
h5write(embedding, file=f, name="/embedding")
h5write(scatterer, file=f, name="/scatterer")
h5write(computation, file=f, name='/computation')

## write attributes
fid <- H5Fopen(f)
# root level
h5writeAttribute("Au prolate spheroid in water", fid, "name")
h5writeAttribute("Computation using SMARTIES, a numerically robust EBCM implementation for spheroids", fid, "description")
h5writeAttribute("gold, spheroid, ebcm, passive, reciprocal, czinfinity, mirrorxyz", fid, "keywords")
h5writeAttribute("v1", fid, "storage_format_version")

# wavelength
did <- H5Dopen(fid, "vacuum_wavelength")
h5writeAttribute("nm", did, "unit")
H5Dclose(did)

# embedding
gid <- H5Gopen(fid, "embedding")
h5writeAttribute("H2O, Water", gid, "name")
h5writeAttribute("non-dispersive", gid, "keywords")
H5Gclose(gid)

# material
gid <- H5Gopen(fid, "scatterer/material")
h5writeAttribute("Au, Gold", gid, "name")
h5writeAttribute("dispersive, plasmonic", gid, "keywords")
h5writeAttribute("Au from Raschke et al 10.1103/PhysRevB.86.235147", gid, "reference")
H5Gclose(gid)

# geometry
gid <- H5Gopen(fid, "scatterer/geometry")
h5writeAttribute("nm", gid, "unit")
h5writeAttribute("spheroid", gid, "shape")
h5writeAttribute("homogeneous spheroid with symmetry axis z", gid, "name")
H5Gclose(gid)

# computation
gid <- H5Gopen(fid, "computation")
h5writeAttribute("EBCM, Extended Boundary Condition Method", gid, "method")
h5writeAttribute("Computation using SMARTIES, a numerically robust EBCM implementation for spheroids", gid, "description")
h5writeAttribute(software, gid, "software")
h5writeAttribute("SMARTIES", gid, "name")
H5Gclose(gid)

H5Fclose(fid)

