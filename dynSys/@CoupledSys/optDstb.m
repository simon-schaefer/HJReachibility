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