function obj = init( obj, rtab, ztab, ytab, varargin )
%  INIT - Initialization of Green function plotter.
%
%  Usage for obj = stratified.plot :
%    obj = init( obj, rtab, ztab, ytab, PropertyPairs  )
%    rtab     :  radii for tabulated Green functions
%    ztab     :  z-values for tabulated Green functions    
%    ytab     :  tabulated Green function  
%  PropertyPair
%    fun      :  @real, @imag, @abs, or user-defined function

%  set up parser
p = inputParser;
p.KeepUnmatched = true;
addParameter( p, 'fun', @abs );
%  parse input
parse( p, varargin{ : } );

%  save input
obj.rtab = rtab;
obj.ztab = ztab;
obj.ytab = ytab;
%  Green function element and display function for plotting
names = fieldnames( ytab );
obj.name = names{ 1 };
obj.fun = p.Results.fun;

%  figure
figure
imagesc( rtab, ztab, obj.fun( ytab.( obj.name ).Values ) .' );

set( gca, 'YDir', 'norm' );
xlabel( 'r (nm)' );
ylabel( 'Z (nm)' );

set( gcf, 'UserData', obj );
figname( obj );
contextmenu( obj );
initpaging( obj );
