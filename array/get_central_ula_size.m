function [s, m_v, zero_idx] = get_central_ula_size(diffs)
%GET_CENTRAL_ULA_SIZE Obtains the size of the ULA centered at the origin
%from the difference coarray.
%Syntax:
%   s = GET_CENTRAL_ULA_SIZE(diffs);
%   [s, m_v, zero_idx] = GET_CENTRAL_ULA_SIZE(diffs);
%Input:
%   diffs - A vector of integer differences, or array design.
%Outputs:
%   s - Entire array size = 2*M_v - 1.
%   m_v - The size of the ULA after trimming the negative part.
%   zero_idx - Index of the zero element in diffs.

if isstruct(diffs)
    % convert array design to unique differences
    diffs = unique_differences(diffs.element_indices);
end

w = 0;
zero_idx = find(diffs == 0);
n_diff = length(diffs);
while zero_idx + w < n_diff && zero_idx - w > 1
    w = w + 1;
    if (diffs(zero_idx+w) - diffs(zero_idx+w-1) ~= 1) || (diffs(zero_idx-w+1) - diffs(zero_idx-w)) ~= 1
        w = w - 1;
        break;
    end
end
m_v = w + 1;
s = 2*m_v - 1;

end

