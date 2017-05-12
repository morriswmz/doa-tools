clear();

f_number = @(x, y) x*x;
f_matrix = @(x, y) [x x*x];
f_string = @(x, y) num2str(x);
f_struct = @(x, y) struct('value', x);

n = 16;
r = 20;
params = [(1:n)' zeros(n,1)];

%% No parallel execution. Function returns scalar.
val_expected = repmat(params(:,1).^2, 1, r);
val_actual = run_for_parameters(f_number, params, 'Repeat', r);
assert(isequal(val_actual, val_expected));

%% No parallel execution. Function returns matrix.
val_expected = cell(n, r);
for ii = 1:n
    val_expected(ii,:) = {[params(ii,1) params(ii,1)^2]};
end
val_actual = run_for_parameters(f_matrix, params, 'Repeat', r);
assert(isequal(val_actual, val_expected));

%% No parallel execution. Function returns string.
val_expected = cell(n, r);
for ii = 1:n
    val_expected(ii,:) = {num2str(params(ii,1))};
end
val_actual = run_for_parameters(f_string, params, 'Repeat', r);
assert(isequal(val_actual, val_expected));

%% No parallel execution. Function returns struct.
val_expected = repmat(struct('value', []), n, 1);
for ii = 1:n
    val_expected(ii) = struct('value', params(ii,1));
end
val_expected = repmat(val_expected, 1, r);
val_actual = run_for_parameters(f_struct, params, 'Repeat', r);
assert(isequal(val_actual, val_expected));

%% Parallel level = 0. Function returns scalar.
val_expected = repmat(params(:,1).^2, 1, r);
val_actual = run_for_parameters(f_number, params, 'Repeat', r, 'Parallel', true);
assert(isequal(val_actual, val_expected));

%% Parallel level = 0. Function returns matrix.
val_expected = cell(n, r);
for ii = 1:n
    val_expected(ii,:) = {[params(ii,1) params(ii,1)^2]};
end
val_actual = run_for_parameters(f_matrix, params, 'Repeat', r, 'Parallel', true);
assert(isequal(val_actual, val_expected));

%% Parallel level = 0. Function returns string.
val_expected = cell(n, r);
for ii = 1:n
    val_expected(ii,:) = {num2str(params(ii,1))};
end
val_actual = run_for_parameters(f_string, params, 'Repeat', r, 'Parallel', true);
assert(isequal(val_actual, val_expected));

%% Parallel level = 0. Function returns struct.
val_expected = repmat(struct('value', []), n, 1);
for ii = 1:n
    val_expected(ii) = struct('value', params(ii,1));
end
val_expected = repmat(val_expected, 1, r);
val_actual = run_for_parameters(f_struct, params, 'Repeat', r, 'Parallel', true);
assert(isequal(val_actual, val_expected));

%% Parallel level = 1. Function returns scalar.
val_expected = repmat(params(:,1).^2, 1, r);
val_actual = run_for_parameters(f_number, params, 'Repeat', r, 'Parallel', true, 'ParallelLevel', 1);
assert(isequal(val_actual, val_expected));

%% Parallel level = 1. Function returns matrix.
val_expected = cell(n, r);
for ii = 1:n
    val_expected(ii,:) = {[params(ii,1) params(ii,1)^2]};
end
val_actual = run_for_parameters(f_matrix, params, 'Repeat', r, 'Parallel', true, 'ParallelLevel', 1);
assert(isequal(val_actual, val_expected));

%% Parallel level = 1. Function returns string.
val_expected = cell(n, r);
for ii = 1:n
    val_expected(ii,:) = {num2str(params(ii,1))};
end
val_actual = run_for_parameters(f_string, params, 'Repeat', r, 'Parallel', true, 'ParallelLevel', 1);
assert(isequal(val_actual, val_expected));

%% Parallel level = 1. Function returns struct.
val_expected = repmat(struct('value', []), n, 1);
for ii = 1:n
    val_expected(ii) = struct('value', params(ii,1));
end
val_expected = repmat(val_expected, 1, r);
val_actual = run_for_parameters(f_struct, params, 'Repeat', r, 'Parallel', true, 'ParallelLevel', 1);
assert(isequal(val_actual, val_expected));

%% Progress bar demonstration
fprintf('Expecting total progress = %d\n', n*r);
run_for_parameters(f_number, params, 'Repeat', r, 'ShowProgress', true);
fprintf('Expecting total progress = %d\n', n);
run_for_parameters(f_number, params, 'Repeat', r, 'ShowProgress', true, 'Parallel', true);
fprintf('Expecting total progress = %d\n', n*r);
run_for_parameters(f_number, params, 'Repeat', r, 'ShowProgress', true, 'Parallel', true, 'ParallelLevel', 1);