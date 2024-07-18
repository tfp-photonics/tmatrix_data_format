## mockup data
using Pkg, HDF5

# possibly multiple wavelengths
wavelength = collect(400:50:800)
Nl = length(wavelength)
Lmax = 3
qmax = 2*(Lmax*(Lmax+1)+Lmax)

# dummy 30x30 matrix values for each wavelength
# note the row-major ordering
tdata = transpose(reshape(collect(1:qmax^2) + collect(1:qmax^2) * 1im, (qmax,qmax)))
tmatrix = zeros(ComplexF64,(Nl,qmax,qmax))
for i=1:Nl
    tmatrix[i,:,:] = tdata
end

print(tmatrix[1,1:3,1:3])

l = zeros(Int64, qmax)
m = zeros(Int64, qmax)
s = Vector{String}(undef,qmax)
let
i=1
for li = 1:Lmax
    for mi = -li:li
        for si = ["electric", "magnetic"]
            l[i] = li
            m[i] = mi
            s[i] = si
            i = i+1
        end
    end
end
end


f = "aj.tmat.h5"
ver = Pkg.Operations.Context().env.manifest
h5ver = string(ver[findfirst(v->v.name == "HDF5", ver)].version)
software = "SMARTIES=1.1, julia=$(VERSION), HDF5.jl=$(h5ver)"

h5open(f, "w") do fid

    fid["vacuum_wavelength"] = wavelength
    attributes(fid["vacuum_wavelength"])["unit"] = "nm"

    # set = permutedims(dset, reverse(1:ndims(dset)))
    # https://juliaio.github.io/HDF5.jl/stable/#Language-interoperability-with-row-and-column-major-order-arrays
    fid["tmatrix"] = permutedims(tmatrix, reverse(1:ndims(tmatrix)))

    modes = create_group(fid, "modes") 
    modes["l"] = l
    modes["m"] = m
    modes["polarization"] = s

    embedding = create_group(fid, "embedding") 
    embedding["relative_permittivity"] = 1.33^2
    embedding["relative_permeability"] = 1.0
    attributes(embedding)["name"] = "H2O, Water"
    attributes(embedding)["keywords"] = "non-dispersive"

    sca_mat = create_group(fid, "scatterer/material") 
    sca_mat["relative_permittivity"] = repeat([-11.4 + 1.181im], Nl)
    sca_mat["relative_permeability"] = 1.0
    attributes(sca_mat)["name"] = "Au, Gold"
    attributes(sca_mat)["keywords"] = "dispersive, plasmonic"
    attributes(sca_mat)["reference"] = "Au from Raschke et al 10.1103/PhysRevB.86.235147"

    sca_geo = create_group(fid, "scatterer/geometry") 
    sca_geo["radiusxy"] = 20.0
    sca_geo["radiusz"] = 40.0
    attributes(sca_geo)["unit"] = "nm"
    attributes(sca_geo)["shape"] = "spheroid"
    attributes(sca_geo)["name"] = "homogeneous spheroid with symmetry axis z"

    mpar = create_group(fid, "computation/method_parameters") 
    mpar["Lmax"] = Lmax
    mpar["Ntheta"] = 100
    script = create_group(fid, "computation/files") 
    script["script"] = read("export_julia.jl", String)
    
    # write root attributes
    attributes(fid)["name"] = "Au prolate spheroid in water"
    attributes(fid)["description"] = "Computation using SMARTIES, a numerically robust EBCM implementation for spheroids"
    attributes(fid)["keywords"] = "gold, spheroid, ebcm, passive, reciprocal, czinfinity, mirrorxyz"    
    attributes(fid)["storage_format_version"] = "v1"
    
    # comp attributes
    attributes(fid["computation"])["method"] = "EBCM, Extended Boundary Condition Method"
    attributes(fid["computation"])["description"] = "Computation using SMARTIES, a numerically robust EBCM implementation for spheroids"
    attributes(fid["computation"])["name"] = "SMARTIES"
    attributes(fid["computation"])["software"] = software

end

