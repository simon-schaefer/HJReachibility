function run_coupled_sys()
 
%% Environment parameters
ndim = 2;

x_min = -20;
x_max = 20;
v_min_robot = -2; 
v_max_robot = 2;
a_max_robot = 2;
v_max_single = 4;
N = [81; 81; 21; 21];

%% Grid (isotropic)
grid_min = [ones(ndim,1) * x_min; ones(ndim,1) * v_min_robot]';
grid_max = [ones(ndim,1) * x_max; ones(ndim,1) * v_max_robot]';

g = createGrid(grid_min, grid_max, N);


%% target set
R = 1;
% data0 = shapeCylinder(grid,ignoreDims,center,radius)
data0 = shapeCylinder(g, ndim+1:2*ndim, zeros(ndim), R);

%% time vectors
t0 = 0;
tMax = 2.0;
dt = 0.4;
tau = t0:dt:tMax;

%% problem parameters

% input and disturbance bounds
aMax = a_max_robot;
dMax = ones(ndim) * v_max_single;

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
output_file = sprintf("../reachability/%dD.mat", ndim);
save(output_file, 'grid_min', 'grid_max', 'N', ...
    'value_function', 'gradients', 'tau');
end