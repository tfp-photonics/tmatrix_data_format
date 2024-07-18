## convert a tmat.h5 file using different units to wavelength in nm
# The example below reads a tmat file generated with Treams, where the frequency information is not stored using the field `vacuum_wavelength` in units of nm, as required in TERMS. We therefore copy the file, read the frequency and its unit, convert those to wavelengths in nm and write these new data to the file

using Pkg, HDF5
pwd()
# conversions between units
include("unit_conversions.jl")

# editing a copy of the original file
f1 = "../reference_data/tetrahedron.tmat.h5"
f2 = "../terms/tetrahedron_nm.tmat.h5"
cp(f1, f2, force=true)

fid = h5open(f2, "r+")
g = fid["/"]
keys(g)

if "frequency" in keys(g)
    f = read_dataset(fid, "frequency")
    unit = attributes(fid["frequency"])["unit"][]
    wavelength = convert_frequency(f, unit)
    delete_object(fid, "frequency")
elseif "angular_frequency" in keys(g)
    ω = read_dataset(fid, "angular_frequency")
    unit = attributes(fid["angular_frequency"])["unit"][]
    wavelength = convert_angular_frequency(ω, unit)
    delete_object(fid, "angular_frequency")
elseif "vacuum_wavelength" in keys(g)
    λ = read_dataset(fid, "vacuum_wavelength")
    unit = attributes(fid["vacuum_wavelength"])["unit"][]
    wavelength = convert_wavelength(λ, unit)
    delete_object(fid, "vacuum_wavelength")
elseif "vacuum_wavenumber" in keys(g)
    ν = read_dataset(fid, "vacuum_wavenumber")
    unit = attributes(fid["vacuum_wavenumber"])["unit"][]
    wavelength = convert_wavenumber(ν, unit)
    delete_object(fid, "vacuum_wavenumber")
elseif "angular_vacuum_wavenumber" in keys(g)
    k = read_dataset(fid, "angular_vacuum_wavenumber")
    unit = attributes(fid["angular_vacuum_wavenumber"])["unit"][]
    wavelength = convert_angular_wavenumber(k, unit)
    delete_object(fid, "angular_vacuum_wavenumber")
else
    message("The file does not appear to have correct frequency information")
end

write_dataset(fid, "vacuum_wavelength", round.(wavelength; digits=3)) # rounding errors
attributes(fid["vacuum_wavelength"])["unit"] = "nm"
read_dataset(fid, "vacuum_wavelength")
close(fid)
