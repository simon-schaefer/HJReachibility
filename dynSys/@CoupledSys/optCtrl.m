function uOpt = optCtrl(obj, ~, ~, deriv, uMode)
% uOpt = optCtrl(obj, t, x, deriv, uMode)

%% Input processing
if nargin < 5
  uMode = 'min';
end

if ~iscell(deriv)
  deriv = num2cell(deriv);
end

uOpt = cell(obj.nu, 1);
Rmax = 1;

%% Optimal control
% The avoid set is assumed to be in the center, i.e. around x_i = 0. 
% Then in order to avoid it at every case try to distance yourself from 
% it as best as possible, i.e. when you are positive accelerate in further 
% positive direction, and vice versa.
% This is also the result of the mathematical formulation:
% d/du \nabla V * f(x, d) = dV3 ux + dV4 uy
% so that we maximize this equation when u has the same sign and maximal
% absolute value as dVi (dVi = ith component of value function gradient). 
if strcmp(uMode, 'max')
    for i = 1:obj.nu
        uOpt{i} = (deriv{i}>=0)*(obj.aRange(2))+(deriv{i}<0)*obj.aRange(1);
    end
    
elseif strcmp(uMode, 'min')
    for i = 1:obj.nu
        uOpt{i} = (deriv{i}>=0)*(obj.aRange(2))+(deriv{i}<0)*obj.aRange(1);
    end 
    
else
  error('Unknown uMode!')
end

end