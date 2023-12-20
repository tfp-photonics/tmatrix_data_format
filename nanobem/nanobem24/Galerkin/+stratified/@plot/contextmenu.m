function contextmenu( varargin )
%  CONTEXTMENU - Add context menu to current figure.

%  context menu
h = uicontextmenu( gcf );
%  callback items for function
uimenu( h, 'Label', 'abs',  'Callback', @( ~, ~ ) fun( @abs  ) );
uimenu( h, 'Label', 'real', 'Callback', @( ~, ~ ) fun( @real ) );
uimenu( h, 'Label', 'imag', 'Callback', @( ~, ~ ) fun( @imag ) );

%  attach menu to figure
set( gcf, 'uicontextmenu', h );


function fun( varargin )
%  SET - Set and apply user-defined function.

%  get stratified.plot object and modify function
obj = get( gcf, 'UserData' );
obj.fun = deal( varargin{ : } );

%  refresh figure
set( get( gca, 'Children' ),  ...
  'CData', obj.fun( obj.ytab.( obj.name ).Values ) .' );
figname( obj );
drawnow;
%  save stratified.plot object in figure
set( gcf, 'UserData', obj );
