function approximate()
%% Data loading
ndim = 2;
data_file = "../reachability/1D.mat";
output_file = "../reachability/1D_app.mat";
up_scaling = 8;
training_steps = 10;

% (Re-)Create grid points from data.
load(data_file, 'grid_min', 'grid_max', 'N', 'value_function');
X = makegrid(grid_min, grid_max, N, ndim);
Y = makevalues(value_function, ndim);
Xt = X; 
Yt = Y;
[n, ~] = size(X);

% Create up-scaled grid points for prediction later on.
N = up_scaling*N;
Xp = makegrid(grid_min, grid_max, N, ndim);

tic

%% LWPR initialization
model = lwpr_init(ndim, 1, 'name', 'value_function_approximation');
model = lwpr_set(model,'init_D', eye(ndim)*25);    
model = lwpr_set(model,'init_alpha', ones(ndim)*250);
model = lwpr_set(model,'diag_only',0);
model = lwpr_set(model,'w_gen', 0.2);
% model = lwpr_set(model,'w_prune', 0.7);   
model = lwpr_set(model,'diag_only', 0);   
model = lwpr_set(model,'meta', 1);
model = lwpr_set(model,'meta_rate', 250);
model = lwpr_set(model,'kernel', 'Gaussian');  
% model = lwpr_set(model,'update_D', 1);

% Transfer model into mex-internal storage
model = lwpr_storage('Store', model);

%% Model training
nmse = zeros(training_steps, 1);
f = waitbar(0, 'Training model ...');
for j = 1:training_steps
   inds = randperm(n);

   mse = 0;
   for i=1:n
       waitbar(i/n, f, sprintf('Training (%d / %d) ==> %d / %d ...', ...
                                j, training_steps, i, n));
	   [model, yp, ~] = lwpr_update(model,X(inds(i),:)',Y(inds(i),:)');         
	   mse = mse + (Y(inds(i),:)-yp).^2;
   end

   nMSE = mse/n/var(Y,1);
   fprintf(1, '#Data=%d #rfs=%d nMSE=%5.3f\n', ...
           lwpr_num_data(model), lwpr_num_rfs(model), nMSE);
   %if exist('fflush') % for Octave output only
   %   fflush(1);
   %end   
   nmse(j) = nMSE;
end
close(f);


%% Model predictions
Yp = zeros(size(Xp, 1), size(Yt, 2));
NXp = length(Xp);
f = waitbar(0, 'Predicting ...');
for i=1:NXp
    waitbar(i/NXp, f, sprintf('Predicting %d / %d ...', i, NXp));
	[yp, ~] = lwpr_predict(model, Xp(i,:)', 0.001);
	Yp(i,1) = yp;
end
close(f);
% [yp, ~] = lwpr_predict(model, Xp', 0.001);
% Yp = yp';

toc

%% Model cleanup
% Transfer model back from mex-internal storage
model = lwpr_storage('GetFree', model);

%% Plotting output data
figure(1);
clf;

% plot the raw noisy data
subplot(2,2,1);
plot3(X(:,1),X(:,2),Y,'*');
title('Noisy data samples');

% plot the fitted surface
axis([-1 1 -1 1 -.5 1.5]);
subplot(2,2,2);
[x,y,z]=makesurf([Xp,Yp],sqrt(length(Xt)));
if ~exist('surfl')
   mesh(x,y,z);
else
   surfl(x,y,z);
end
axis([-1 1 -1 1 -.5 1.5]);
title(sprintf('The fitted function: nMSE=%5.3f',nmse));

% plot the true surface
subplot(2,2,3);
[x,y,z]=makesurf([Xt,Yt],sqrt(length(Xt)));
if ~exist('surfl')
   mesh(x,y,z);
else
   surfl(x,y,z);
end
axis([-1 1 -1 1 -.5 1.5]);
title('The true function');

% plot the local models
subplot(2,2,4);
for i=1:length(model.sub(1).rfs)
    % D = R*model.sub(1).rfs(i).D*R';
    % c = R*model.sub(1).rfs(i).c;
    % draw_ellipse(D(1:2,1:2),c(1:2),0.1,model.kernel);
	draw_ellipse(model.sub(1).rfs(i).D, model.sub(1).rfs(i).c, ...
                 0.1,'Gaussian');
	hold on;
end
hold off;
axis('equal');
title('Projected input space view of RFs');
% stitle('Input space view of RFs');

%% Plotting training
figure(2);
plot(log(nmse));
xlabel('Iteration');
ylabel('Log(nMSE)');
title('nMSE');

%% Write results to mat file
value_function = Yp;
value_function_flat = reshape(Yp, 1, []);
save(output_file, 'grid_min', 'grid_max', 'N', ...
    'value_function', 'value_function_flat');


% -------------------------------------------------------------------------
function X = makegrid(grid_min, grid_max, N, ndim)
% [X]=makegrid(grid_min, grid_max, N, ndim) creates an n-dimensional
% meshgrid from equally spaced axes with min-max values given in the 
% arrays grid_min and grid_max and number of points in array N.
grid_data = {};
for d = 1:ndim
    grid_data{d} = linspace(grid_min(d), grid_max(d), N(d));
end
if ndim == 2
    [grid_x, grid_vx] = ndgrid(grid_data{:});
    grid_x = permute(grid_x, [2, 1]);  % ngrid -> meshgrid
    grid_vx = permute(grid_vx, [2, 1]);
    X = [grid_x(:), grid_vx(:)];
elseif ndim == 4
    [grid_x, grid_y, grid_vx, grid_vy] = ndgrid(grid_data{:});
    grid_x = permute(grid_x, [2, 1, 3, 4]);  % ngrid -> meshgrid
    grid_y = permute(grid_y, [2, 1, 3, 4]);
    grid_vx = permute(grid_vx, [2, 1, 3, 4]);
    grid_vy = permute(grid_vy, [2, 1, 3, 4]);
    X = [grid_x(:), grid_y(:), grid_vx(:), grid_vy(:)]; 
else
    error('Undefined grid creation for dimension !');
end


function Y = makevalues(values, ndim)
if ndim == 2
    Y = reshape(values(:, :, end), [], 1);
elseif ndim == 4
    Y = reshape(values(:, :, :, :, end), [], 1);
else
    error('Undefined values creation for dimension !');
end

% -------------------------------------------------------------------------
function [X,Y,Z] = makesurf(data,nx)
% [X,Y,Z]=makesurf(data,nx) converts the 3D data file data into
% three matices as need by surf(). nx tells how long the row of the
% output matrices are

[m, ~] = size(data);

n=0;
X = zeros(nx, m);
Y = zeros(nx, m);
Z = zeros(nx, m);
for i=1:nx:m
	n = n+1;
	X(:,n) = data(i:i+nx-1,1);
	Y(:,n) = data(i:i+nx-1,2);
	Z(:,n) = data(i:i+nx-1,3);
end


% -------------------------------------------------------------------------
function []=draw_ellipse(M,C,w,kernel)
% function draw ellipse draws the ellipse corresponding to the
% eigenvalues of M at the location c.

[V,E] = eig(M);
d1 = E(1,1);
d2 = E(2,2);

steps = 50;
switch kernel
case 'Gaussian'
	start = sqrt(-2*log(w)/d1);
case 'BiSquare'
	start = sqrt(4*(1-sqrt(w))/d1);
end

Xp = zeros(steps, 1);
Yp = zeros(steps, 1);
for i=0:steps
	Xp(i+1,1) = -start + i*(2*start)/steps;
	
    switch kernel
	case 'Gaussian'
		arg = (-2*log(w)-Xp(i+1,1)^2*d1)/d2;
	case 'BiSquare'
		arg = (2*(1-sqrt(w))-Xp(i+1,1)^2*d1)/d2;
    end
    
    if (arg < 0)
		arg = 0; 
    end % should be numerical error
	Yp(i+1,1) = sqrt(arg);
end

for i=1:steps+1
	Xp(steps+1+i,1) = Xp(steps+1-i+1,1);
	Yp(steps+1+i,1) = -Yp(steps+1-i+1,1);
end

% transform the rf

M = [Xp,Yp]*V(1:2,1:2)';

Xp = M(:,1) + C(1);
Yp = M(:,2) + C(2);

plot(C(1),C(2),'ro',Xp,Yp,'c');