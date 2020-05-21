function Y = makeValues(values, ndim)
if ndim == 2
    Y = reshape(values(:, :, end), [], 1);
elseif ndim == 4
    Y = reshape(values(:, :, :, :, end), [], 1);
else
    error('Undefined values creation for dimension !');
end
