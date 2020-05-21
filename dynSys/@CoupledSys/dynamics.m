function dx = dynamics(obj, ~, x, u, d)
% Dynamics of the Coupled System (exampled ndim = 2)
      %    \dot{x}_1 = vx - d1
      %    \dot{x}_2 = vy - d2
      %    \dot{x}_3 = ax
      %    \dot{x}_4 = ay
%   Control: u = ax, ay;

if nargin < 5
  d = zeros(obj.ndim, 1);
end

if iscell(x)
  dx = cell(length(obj.dims), 1);
  for i = 1:length(obj.dims)
    dx{i} = dynamics_cell_helper(obj, x, u, d, obj.dims(i));
  end
else
  dx = zeros(obj.nx, 1);
  
  dx(1:ndim) = x(ndim+1:end) - d;
  dx(ndim+1:end) = u;
end
end

function dx = dynamics_cell_helper(obj, x, u, d, dim)

if dim <= obj.ndim
    dx = x{dim + obj.ndim} - d{dim};
else
    dx = u{dim - obj.ndim};
end

end