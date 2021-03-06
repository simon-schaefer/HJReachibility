function run_coupled_sys()
 
%% Environment parameters
ndim = 2;
N = [41; 41; 21; 21];

x_min = -30;
x_max = 30;
v_min_robot = -2.0; 
v_max_robot = 2.0;
a_max_robot = 2;
v_max_ped = 2.5;

t0 = 0;
tMax = 6.0;
dt = 0.4;

output_file = sprintf("../reachability/%dD.mat", ndim);

%% Grid (isotropic)
grid_min = [ones(ndim,1) * x_min; ones(ndim,1) * v_min_robot]';
grid_max = [ones(ndim,1) * x_max; ones(ndim,1) * v_max_robot]';

g = createGrid(grid_min, grid_max, N);

%% target set
R = 1;
% data0 = shapeCylinder(grid,ignoreDims,center,radius)
data0 = shapeCylinder(g, ndim+1:2*ndim, zeros(ndim), R);

%% time vectors
tau = t0:dt:tMax;

%% problem parameters

% input and disturbance bounds
aMax = a_max_robot;
dMax = ones(ndim) * v_max_ped;

% control trying to min or max value function?
uMode = 'max';
dMode = 'min';
minWidth = 'minVOverTime';

%% Pack problem parameters

% Define dynamic system
dSys = CoupledSys(zeros(2*ndim,1), ndim, aMax, dMax);

% Put grid and dynamic systems into schemeData
schemeData.grid = g;
schemeData.dynSys = dSys;
schemeData.accuracy = 'high'; %set accuracy
schemeData.uMode = uMode;
schemeData.dMode = dMode;

%% Compute value function

%HJIextraArgs.visualize = true; %show plot
HJIextraArgs.visualize.valueSet = 1;
HJIextraArgs.visualize.initialValueSet = 1;
HJIextraArgs.visualize.figNum = 1; %set figure number
HJIextraArgs.visualize.deleteLastPlot = true; %delete previous plot

%[data, tau, extraOuts] = ...c
% HJIPDE_solve(data0, tau, schemeData, minWith, extraArgs)
[value_function, tau, ~] = ...
  HJIPDE_solve(data0, tau, schemeData, minWidth, HJIextraArgs);
[gradients, ~, ~] = computeGradients(g, value_function);

%% Write results to mat file
params = struct('v_max_ped', v_max_ped, ...
                'v_max_robot', v_max_robot, ...
                'a_max_robot', a_max_robot);
            
save(output_file, 'grid_min', 'grid_max', 'N', ...
    'value_function', 'gradients', 'tau', 'params');
end