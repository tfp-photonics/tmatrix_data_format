## read long data
# The example below reads T-matrix data from a plain text file output by either SMARTIES (function 
# `exportTmatrix()`) or TERMS (keyword `DumpCollectiveTmatrix`). 

# ```
#   s sp l lp  m mp            Tr            Ti
# 1 1  1 1  1 -1 -1 -6.049214e-05 -4.266526e-04
# 2 1  1 1  1  0  0 -3.331557e-05 -3.932179e-04
# 3 1  1 1  1  1  1 -6.049214e-05 -4.266526e-04
# 4 1  1 1  3 -1 -1 -2.374705e-07 -1.995117e-06
# 5 1  1 1  3  0  0 -1.110299e-07 -1.278537e-06
# 6 1  1 1  3  1  1 -2.374705e-07 -1.995117e-06
# ...
# ```

using Pkg, HDF5, DelimitedFiles, DataFrames, CSV
pwd()

# read all the comments to extract column names, wavelength and epsilon
f = readlines(open("smarties.txt"))
breaks = findall(occursin.(r"^#",f))
f[breaks]
reps = diff([breaks[2:end]; length(f)+1])
# meta = filter(line -> occursin(r"^#",line),f)
meta = f[breaks]
colnames = string.(split(meta[1], ' '))[2:9]
m1 = match.(r"lambda= ([0-9]+)",meta[2:end])
m2 = match.(r"(?<=epsIn= ).*$",meta[2:end])
wavelength = map(x -> parse(Float64,x.captures[1]), m1)
epsilon = map(x -> parse(Complex{Float64},x.match), m2)

Nl = length(wavelength)

# read the T-matrix data
d = CSV.read("smarties.txt", DataFrame, delim=' ', ignorerepeated=true, comment="#", header=colnames)
# split by wavelength
wv = reduce(vcat, [repeat([wavelength[i]], n-1) for (i,n) in enumerate(reps)])
d.wavelength = wv
dg = groupby(d, :wavelength, sort=false)

Lmax = maximum(d.n)
qmax = 2*(Lmax*(Lmax+1)+Lmax)

# indices u, u' to place the T-matrix elements
# in alternating electric/magnetic order

tmat_indexing = function(l,m,s,lmax)
    p = l * (l + 1) + m
    pmax = lmax*(lmax+1)+lmax
    q = (s-1) * pmax + p
    return 2*(p-1) + (3 - s) # magnetic/electric -> electric/magnetic
end

tmatrix = zeros(ComplexF64,(Nl,qmax,qmax))
for i in 1:Nl

u = tmat_indexing.(dg[i].n, dg[i].m, dg[i].s,Lmax)
up = tmat_indexing.(dg[i].np, dg[i].mp, dg[i].sp,Lmax)
    tmatrix[CartesianIndex.(i, u, up)] .= dg[i].Tr + 1im*dg[i].Ti
end

l = zeros(Int64, qmax)
m = zeros(Int64, qmax)
s = Vector{String}(undef,qmax)
let
i=1
for li = 1:Lmax
    for mi = -li:li
        for si = ["electric","magnetic"]
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
    sca_mat["relative_permittivity"] = epsilon
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
    script["script"] = read("convert_julia.jl", String)
    
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

