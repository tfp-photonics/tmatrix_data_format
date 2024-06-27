classdef laguerregauss
  %  Laguerre-Gauss beam in paraxial approximation.
  
  properties
    mat       %  material function
    foc       %  effective focal length
    w0        %  input beam waist
    pol       %  polarization before focusing
    nLG       %  radial index
    mLG       %  azimuthal index
    pow = 1   %  laser power (W)
    M2  = 1   %  M-squared value
  end  
  
  properties (Hidden, Dependent)
    fun       %  Laguerre function
  end  
  
  methods
    function obj = laguerregauss( varargin )
      %  Initialize Laguerre-Gauss beam.
      %
      %  Usage :
      %    obj = LGfield( mat, PropertyPairs )
      %  Input
      %    mat    :  material properties of embedding medium
      %  PropertyName
      %    foc    :  effective focal length
      %    w0     :  input beam waits
      %    pol    :  polarization before focusing
      %    nLG    :  radial index
      %    mLG    :  azimuthal index
      %    pow    :  laser power (W)
      %    M2     :  M-squared value     
      if ~isempty( varargin ),  obj = init( obj, varargin{ : } );  end
    end
    
    function varargout = subsref( obj, s )
      %  Evaluate Lagurre-Gauss object.
      switch s( 1 ).type
        case '()'
          [ varargout{ 1 : nargout } ] = eval( obj, s( 1 ).subs{ : } );
        otherwise
          [ varargout{ 1 : nargout } ] = builtin( 'subsref', obj, s );
      end    
    end
    
    function w1 = waist( obj, k0 )
      %  WAIST - Beam waist in focus.
      w1 = 2 * obj.foc * obj.M2 / ( k0 * obj.w0 );
    end
    
    function wmax = wmax( obj, k0, z )
      %  WMAX - Waist at maximum.
      w1 = waist( obj, k0 );
      zr = 0.5 * obj.mat.k( k0 ) * w1 ^ 2;
      wmax = sqrt( max( obj.mLG / 2, 0.5 ) ) * w1 * sqrt( 1 + ( z / zr ) .^ 2 );
    end
    
    function fun = get.fun( obj )
      %  FUN - Laguerre functio, unnormalized.
      fun = @( x ) x .^ obj.mLG .*  ...
             laguerre( obj.nLG, obj.mLG, x .^ 2 ) .* exp( - 0.5 * x .^ 2 );
    end    
  end
end
