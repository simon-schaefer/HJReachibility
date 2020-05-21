function approximate_interpn()
%% Data loading
data_file = "../reachability/2D.mat";
output_file = "../reachability/2D_app.mat";
k = 2;
ndims = 2;  % equivalent to ndims in run_coupled_sys()
method = 'spline';

%% Use interpolation for up-scaling for value function and gradients.
load(data_file, 'N', 'value_function', 'gradients');
time_steps = size(value_function, 2*ndims + 1);

% Value function.
N_with_time  = [N * k - 1; ones(1, 1) * time_steps];
value_function_app = zeros(N_with_time');
size(value_function_app)
for i = 1:time_steps
    V = value_function(:, :, :, :, i);
    value_function_app(:, :, :, :, i) = interpn(V, method);
end

% Gradients.
N_with_time  = [ones(1, 1) * ndims * 2; N_with_time];
gradients_app = zeros(N_with_time');
size(gradients_app)
for i = 1:time_steps
    for j = 1:ndims
        grad = gradients{j}(:, :, :, :, i);
        gradients_app(j, :, :, :, :, i) = interpn(grad, method);
    end
end

%% Reduce precision (reducing loading time).
disp('Reducing size (float -> single precision) ...');
value_function = single(value_function_app); 
gradients = single(gradients_app);

%% Write results to mat file (version 7.3 for mat-files > 2GB).
disp('Storing results ...');
load(data_file, 'grid_min', 'grid_max', 'tau');
save(output_file, 'grid_min', 'grid_max', 'N', ...
    'value_function', 'gradients', 'tau', '-v7.3');
