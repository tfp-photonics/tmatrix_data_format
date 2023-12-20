function k0 = to_wavenumber( val, name, unit )
%  TO_WAVENUMBER - Convert input to wavenumber.
%
%  Usage :
%    k0 = multipole.to_wavenumber( val, name, unit )
%  Input
%    val    :  value 
%    name   :  "frequency", "vacuum_wavelength", "vacuum_wavenumber"
%    unit   :  units
%  Output
%    k0     :  wavenumber of light in vacuum

%  length units
key = [ "nm", "um", "mm", "cm", "m" ];
value = [ 1e-9, 1e-6, 1e-3, 1e-2, 1 ];
%  length maps
M1 = containers.Map( key, value );
M2 = containers.Map( arrayfun( @( x ) x + "^{-1}", key, 'uniform', 1 ), 1 ./ value );
%  frequencies
key = [ "Hz", "kHz", "MHz", "GHz", "THz" ];
value = [ 1, 1e3, 1e6, 1e9, 1e12 ];
%  frequency map
M3 = containers.Map( key, value );

%  speed of light
clight = 299792458;
%  convert to wavenumber in SI
switch name
  case "frequency"
    k0 = 2 * pi * val * M3( unit ) / clight;
  case "angular_frequency"
    k0 = val * M3( unit ) / clight;
  case "vacuum_wavelength"
    k0 = 2 * pi ./ ( val * M1( unit ) );
  case "vacuum_wavenumber"
    k0 = 2 * pi * val * M2( unit );
  case "angular_vacuum_wavenumber"
    k0 = val * M2( unit );
end
%  convert to nm^{-1}
k0 = 1e-9 * k0;
