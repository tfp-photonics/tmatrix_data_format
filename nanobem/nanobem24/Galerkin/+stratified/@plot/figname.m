function figname( obj )
%  FIGNAME - Set name of figure.
%
%  Usage for obj = stratified.plot :
%    figname( obj )

set( gcf, 'Name', strjoin( { obj.name, ' | ', func2str( obj.fun ) } ) );
