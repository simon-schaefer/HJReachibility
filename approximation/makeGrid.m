function X = makeGrid(grid_min, grid_max, N, ndim)
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
