function info = h5info( varargin )
%  H5INFO - Addtional information for storage of T-matrices.
%
%  Usage :
%    info = tmatrix.h5info( PropertyPairs )
%  PropertyName
%    tau              :  discretized particle boundary
%    name             :  concisive description of T-matrix
%    description      :  description, e.g. optimization goal
%    matgroupname     :  group name for materials
%    matname          :  material name
%    matdescription   :  material description
%    files            :  save Matlab files 
%  Output
%    info             :  structure with additional information

%  set up parser
p = inputParser;
p.KeepUnmatched = true;
addOptional( p, 'tau', [] );
addParameter( p, 'name', "" );
addParameter( p, 'description', "" );
addParameter( p, 'matgroupname', "" );
addParameter( p, 'matname', [] );
addParameter( p, 'matdescription', [] );
addParameter( p, 'files', [] );
%  parse input
parse( p, varargin{ : } );

for name = convertCharsToStrings( fieldnames( p.Results ) ) .'
  info.( name ) = p.Results.( name );
end
