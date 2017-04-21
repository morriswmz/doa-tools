function F = coarray_selection_matrix_1d(design, exclude_negative_part, C)
%COARRAY_SELECTION_MATRIX_1D Construct the coarray selection matrix for
%extracting the virtual measurement vector z.
%Syntax:
%   F = COARRAY_SELECTION_MATRIX_1D(design[, exclude_negative_part, C]);
%Inputs:
%   design - Array design.
%   exclude_negative_part - (Optional) if set to true, negative part of the
%                           central ULA will be excluded, and the resulting
%                           F will be a M_v x M^2 matrix.
%   C - (Optional) a M x M matrix that applies custom weights over each
%       element in the M x M sample covariance matrix R.
%Output:
%   F - The coarray selection matrix.
if nargin < 2
    exclude_negative_part = false;
end
if nargin < 3
    C = [];
else
    if any(size(C) ~= [design.element_count design.element_count])
        error('The dimension of C does not matrix that of R.');
    end
end

[diffs, weights, indices] = weight_function_1d(design);
[~, m_v, zero_idx] = get_central_ula_size(diffs);
if exclude_negative_part
    unif_idx_range = zero_idx:zero_idx+m_v-1;
else
    unif_idx_range = zero_idx-m_v+1:zero_idx+m_v-1;
end

F = zeros(length(unif_idx_range), design.element_count^2);
if isempty(C)
    % redundancy averaging
    for i = 1:length(unif_idx_range)
       cur_weight = weights(unif_idx_range(i));
       F(i, indices{unif_idx_range(i)}) = 1/cur_weight;
    end
else
    % custom weight over each element in R
    c = C(:);
    c(c == 0) = eps*10;
    for i = 1:length(unif_idx_range)
       cur_indices = indices{unif_idx_range(i)};
       F(i,cur_indices) = c(cur_indices);
    end
    % normalize each row
    F = bsxfun(@times, F, 1./sum(F,2));
end

end

