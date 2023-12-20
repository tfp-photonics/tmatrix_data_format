function rel = relation( i1, i2 )
%  RELATION - Relation between observation and source media.
%
%  Usage :
%    rel = relation( i1, i2 )
%  Input
%    i1     :  medium index of observation point
%    i2     :  medium index of source point
%  Output
%    rel    :  'same', 'above', 'below'

if i1 == i2
  rel = 'same';
elseif i1 > i2
  rel = 'above';
else
  rel = 'below';
end
