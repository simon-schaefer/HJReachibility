function dOpt = optDstb(obj, ~, x, ~, dMode)
% dOpt = optCtrl(obj, t, y, deriv, dMode)

%% Input processing
if nargin < 5
  dMode = 'max';
end

if ~iscell(x)
  x = num2cell(x);
end

dOpt = cell(obj.nd, 1);

%% Optimal control
if strcmp(dMode, 'max')
    for i = 1:obj.nd
        dOpt{i} = (x{i}>=0)*obj.dRange{2}(i)+(x{i}<0)*(obj.dRange{1}(i));
    end
    
elseif strcmp(dMode, 'min')
    for i = 1:obj.nd
        dOpt{i} = (x{i}>=0)*(obj.dRange{1}(i))+(x{i}<0)*obj.dRange{2}(i);
    end
  
else
  error('Unknown dMode!')
end

end