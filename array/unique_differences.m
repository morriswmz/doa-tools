function diffs = unique_differences(pos)
%UNIQUE_DIFFERENCES Computes the unique position differences.
%Syntax:
%   diffs = UNIQUE_DIFFERENCES([1 2 4 8 9]);
%Input:
%   pos - A matrix where each row represents a position coordinate.
%Output:
%   diffs - A matrix containing all the unique position differences, where
%           each row represents a position coordinate.
[n, ~] = size(pos);
if isvector(pos) && n == 1
    pos = pos(:);
    n = length(pos);
end
p1 = repmat(pos, n, 1);
p2 = kron(pos, ones(n, 1));
diffs = unique(p1 - p2, 'rows');
end

