
c  = 299792458 #m/s

units = Dict(
"yHz"  => 1e-24 , "ym"  => 1e-24 ,         
"zHz"  => 1e-21 , "zm"  => 1e-21 ,        
"aHz"  => 1e-18 , "am"  => 1e-18 ,         
"fHz"  => 1e-15 , "fm"  => 1e-15 ,        
"pHz"  => 1e-12 , "pm"  => 1e-12 ,         
"nHz"  => 1e-9  , "nm"   => 1e-9  ,       
"uHz"  => 1e-6  , "um"  => 1e-6  ,        
"mHz"  => 1e-3  , "mm"  => 1e-3  ,        
"cHz"  => 1e-2  , "cm"  => 1e-2  ,        
"dHz"  => 1e-1  , "dm"  => 1e-1  ,        
"Hz"   => 1     , "m"   => 1     ,   
"daHz" => 1e1   , "dam" => 1e1   ,        
"hHz"  => 1e2   , "hm"  => 1e2   ,       
"kHz"  => 1e3   , "km"  => 1e3   ,       
"MHz"  => 1e6   , "Mm"  => 1e6   ,        
"GHz"  => 1e9   , "Gm"  => 1e9   ,        
"THz"  => 1e12  , "Tm"  => 1e12  ,         
"PHz"  => 1e15  , "Pm"  => 1e15  ,         
"EHz"  => 1e18  , "Em"  => 1e18  ,         
"ZHz"  => 1e21  , "Zm"  => 1e21  ,         
"YHz"  => 1e24  , "Ym"  => 1e24  ,
)

# standardise s^-1 into Hz
factor = units[replace("Zs^{-1}",  "s^{-1}" => "Hz")]


function convert_frequency(f, unit="THz") 
  cleanedunit = replace(unit,  "s^{-1}" => "Hz")
  s = units[cleanedunit]
  return c ./ (f * s) * 1e9 # in nm
end

# convert_frequency(1.0)

function convert_angular_frequency(ω, unit="THz") 
    cleanedunit = replace(unit,  "s^{-1}" => "Hz")
    s = units[cleanedunit]
    return 2π * c ./(ω * s) * 1e9 # in nm
end

# convert_angular_frequency(2.3, "THz") 

function convert_wavenumber(ν, unit="cm^{-1}") 
    cleanedunit = replace(unit,  "µm" => "um")
    cleanedunit = replace(unit,  "^{-1}" => "")# strip inverse
    s = units[cleanedunit]
    return s ./ ν * 1e9 # in nm
end


# convert_wavenumber(300)

function convert_angular_wavenumber(k, unit="cm^{-1}") 
    return 2π * convert_wavenumber(k, unit) 
end

# convert_angular_wavenumber(300)

function convert_wavelength(λ, unit="nm") 
    cleanedunit = replace(unit,  "µm" => "um")
    s = units[cleanedunit]
    return s .* λ * 1e9 # in nm
end
