function obj = init( obj, varargin )
%  Initialize multipole solution.
%
%  Usage for obj = multipole.solution :
%    obj = init( obj, mat, k0, a, b, ai, bi )
%  Input
%    mat    :  material properties of embedding medium
%    k0     :  wavenumber of light in vacuum
%    a        :  TM scattering coefficients
%    b        :  TE scattering coefficinets
%    ai       :  TM incoming coefficients
%    bi       :  TE incoming coefficinets   
      
%  set up parser
p = inputParser;
p.KeepUnmatched = true;
addOptional( p, 'mat', [] );
addOptional( p, 'k0', [] );
addOptional( p, 'a',  [] );
addOptional( p, 'b',  [] );
addOptional( p, 'ai', [] );
addOptional( p, 'bi', [] );
%  parse input
parse( p, varargin{ : } );
      
for name = convertCharsToStrings( fieldnames( p.Results ) ) .'
  obj.( name ) = p.Results.( name );
end
