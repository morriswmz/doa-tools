function [diffs, w, index_map] = weight_function_1d(design)
%WEIGHT_FUNCTION_1D Compute the weight function of 1D array.
% Syntax: [diffs, w, idx] = WEIGHT_FUNCTION_1D(design);
% Inputs:
%   design - 1D array design or an array of indices.
% Outputs:
%   diffs - A column vector, consists of all possible differences
%           (of indices, integers).
%   w - A column vector, consists of weights corresponding to diff vector.
%   index_map - A cell array matching the difference vector. index_map{i}
%               contains all indices in the vectorized R corresponding
%               to the difference diff(i).
if isstruct(design)
    if design.dim > 1
        error('1D array expected.');
    end
    arr_indices = design.element_indices(:);
elseif isvector(design)
    arr_indices = design(:);
else
    error('Unexpected design format.');
end
diffs = [];
w = [];
index_map = {};
if isempty(arr_indices)
    return;
end
n = length(arr_indices);
P = repmat(arr_indices', n, 1);
D = P' - P; % D_ij = d_i - d_j
[all_diffs, idx] = sort(D(:));
cur_diff = all_diffs(1);
start_idx = 1;
jj = 1;
for ii = 2:(n*n)
    if all_diffs(ii) ~= cur_diff
        diffs(jj) = cur_diff;
        index_map{jj} = idx(start_idx:ii-1);
        w(jj) = ii - start_idx;
        start_idx = ii;
        cur_diff = all_diffs(ii);
        jj = jj + 1;
    end
end
% add the last one
diffs(jj) = cur_diff;
index_map{jj} = idx(start_idx:end);
w(jj) = n*n - start_idx + 1;
end

