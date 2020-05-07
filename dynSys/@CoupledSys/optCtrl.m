function uOpt = optCtrl(obj, ~, x, ~, uMode)
% uOpt = optCtrl(obj, t, x, deriv, uMode)

%% Input processing
if nargin < 5
  uMode = 'min';
end

if ~iscell(x)
  x = num2cell(x);
end

uOpt = cell(obj.nu, 1);
Rmax = 1;

%% Optimal control
% The avoid set is assumed to be in the center, i.e. around x_i = 0. 
% Then in order to avoid it at every case try to distance yourself from 
% it as best as possible, i.e. when you are positive accelerate in further 
% positive direction, and vice versa.
if strcmp(uMode, 'max')
    for i = 1:obj.nu
        uOpt{i} = (x{i}>=0)*(obj.aRange(2))+(x{i}<0)*obj.aRange(1);
    end
    
elseif strcmp(uMode, 'min')
    for i = 1:obj.nu
        uOpt{i} = (x{i}>=0)*(obj.aRange(2))+(x{i}<0)*obj.aRange(1);
    end 
    
else
  error('Unknown uMode!')
end

end