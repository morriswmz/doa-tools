function [ind, mask] = diagonal_indices(n)
%DIAGONAL_INDICES Returns the indices of the diagonal elements of a sqaure
%matrix such that S(ind) = diag(S).
%Syntax:
%   ind = DIAGONAL_INDICES(n);
%   [ind, mask] = DIAGONAL_INDICES(n);
%Input:
%   n - Dimension.
%Outputs:
%   ind - Indices.
%   mask - A boolean vector such that S(mask) = diag(S).
i = (1:n)';
ind = i.*(i+1)./2;
if nargout == 2
    mask = false(n*n, 1);
    mask(ind) = true;
end
end

