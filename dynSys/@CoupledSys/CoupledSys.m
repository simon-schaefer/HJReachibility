classdef CoupledSys < DynSys
  properties  
    % Number of dimensions in position space (1D, 2D, 3D, ...)
    ndim
      
    % Acceleartion
    aRange
      
    % Disturbance
    dRange
    
    % Dimensions that are active
    dims
  end
  
  methods
    function obj = CoupledSys(x, ndim, aRange, dRange, dims)
      % obj = CoupledSys2D(x, wMax, speed, dMax, dims)
      %     Two-Agent coupled system 
      %
      % Dynamics (exampled ndim = 2):
      %    \dot{x}_1 = vx - d1
      %    \dot{x}_2 = vy - d2
      %    \dot{x}_3 = ax
      %    \dot{x}_4 = ay
      %         u \in [-aMax, aMax]
      %         d \in [-dMax, dMax]
      %
      % Inputs:
      %   x      - state: [xpos; ypos]
      %   ndim   - number of dimensions (usually 1, 2)
      %   aRange - maximum acceleration
      %   dMax   - disturbance  nds
      %
      % Output:
      %   obj       - a CoupledSys object
      
      if numel(x) ~= obj.nx
        error('Initial state does not have right dimension!');
      end
      
      if ~iscolumn(x)
        x = x';
      end
      
      if nargin < 3
        aRange = [-1 1];
      end
       
      if nargin < 4
        dRange = {zeros(ndim); zeros(ndim)};
      end
      
      if nargin < 5
        dims = 1:2*ndim;
      end
 
      if numel(aRange) < 2
        aRange = [-aRange; aRange];
      end
      
      if ~iscell(dRange)
        dRange = {-dRange, dRange};
      end
      
      % Basic vehicle properties
      obj.pdim = 1:ndim; % Position dimensions
      obj.vdim = ndim+1:2*ndim; % Velocity dimensions
      obj.ndim = ndim;
      obj.nx = length(dims);
      obj.nu = ndim;
      obj.nd = ndim;
      
      obj.x = x;
      obj.xhist = obj.x;
      
      obj.dRange = dRange;
      obj.aRange = aRange;
      obj.dims = dims;
    end
    
  end % end methods
end % end classdef
