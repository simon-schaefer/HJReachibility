function dOpt = optDstb(obj, ~, ~, deriv, dMode)
% dOpt = optCtrl(obj, t, y, deriv, dMode)

%% Input processing
if nargin < 5
  dMode = 'max';
end

if ~iscell(deriv)
  deriv = num2cell(deriv);
end

dOpt = cell(obj.nd, 1);

%% Optimal control
% The result of the mathematical formulation:
% d/dd \nabla V * f(x, d) = dV1 (-d1) + dV2 (-d2)
% so that we maximize this equation when d has the opposite sign and 
% maximal absolute value as dVi (dVi = ith component of value function 
% gradient). 
if strcmp(dMode, 'max')
    for i = 1:obj.nd
        dOpt{i} = (deriv{i}>=0)*obj.dRange{1}(i)+(deriv{i}<0)*(obj.dRange{2}(i));
    end
    
elseif strcmp(dMode, 'min')
    for i = 1:obj.nd
        dOpt{i} = (deriv{i}>=0)*(obj.dRange{2}(i))+(deriv{i}<0)*obj.dRange{1}(i);
    end
  
else
  error('Unknown dMode!')
end

end