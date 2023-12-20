function obj = initpaging( obj )
%  INITPAGING - Set up figure for paging.

%  green arrow to right
icon = fullfile( matlabroot, '/toolbox/matlab/icons/greenarrowicon.gif' );
[ cdata, map ] = imread( icon );
% convert white pixels into a transparent background
map( sum( map, 2 ) == 3 ) = NaN;
% Convert to RGB space
icon2 = ind2rgb( cdata, map );
icon1 = icon2( :, end : -1 : 1, : );

%  get toolbar icons
tbh = findall( gcf, 'Type', 'uitoolbar' );
%  add left arrow to 
uipushtool( tbh, 'CData', icon1, 'Separator','on',  ...
        'HandleVisibility','off', 'TooltipString', 'Previous field',  ...
        'ClickedCallback', @( ~, ~ ) pagedn );
%  add right arrow
uipushtool( tbh, 'CData', icon2,  ...
        'HandleVisibility','off', 'TooltipString', 'Next field',  ...
        'ClickedCallback', @( ~, ~ ) pageup );
      
      
function [] = pagedn
%  PAGEDN - Page down Green functions.

%  get stratified.plot object
obj = get( gcf, 'UserData' );
%  Green function index
names = fieldnames( obj.ytab );
it = find( strcmp( obj.name, names ) );
%  first array element ?
if it == 1,  return;  end
   
%  refresh figure
obj.name = names{ it - 1 };
set( get( gca, 'Children' ),  ...
  'CData', obj.fun( obj.ytab.( obj.name ).Values ) .' );
figname( obj );
drawnow;
%  save stratified.plot object in figure
set( gcf, 'UserData', obj );

      
function [] = pageup
%  PAGEUP - Page up Green functions.

%  get stratified.plot object
obj = get( gcf, 'UserData' );
%  Green function index
names = fieldnames( obj.ytab );
it = find( strcmp( obj.name, names ) );
%  first array element ?
if it == numel( names ),  return;  end
   
%  refresh figure
obj.name = names{ it + 1 };
set( get( gca, 'Children' ),  ...
  'CData', obj.fun( obj.ytab.( obj.name ).Values ) .' );
figname( obj );
drawnow;
%  save stratified.plot object in figure
set( gcf, 'UserData', obj );
