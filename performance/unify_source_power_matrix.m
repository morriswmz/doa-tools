function P = unify_source_power_matrix(p, n)
%UNIFY_SOURCE_POWER_MATRIX A utility function that unifies scalar,
%vector, or matrix inputs of the source power parameter into a matrix.
%Syntax:
%   P = UNIFY_SOURCE_POWER_MATRIX(p, n);
%Inputs:
%   p - A scalar, vector or matrix.
%       If p is a scalar, returns P = pI.
%       If p is a vector, checks its length and returns P = diag(p).
%       If p is a matrix, checks its dimension and returns P.
%   n - Dimension of the final matrix.
%Outputs:
%   P - Unified source power matrix.
if isscalar(p)
    P = p * eye(n);
elseif isvector(p)
    if length(p) ~= n
        error('The length of the source power vector does not match the source number.');
    end
    P = diag(p);
else
    if any(size(p) ~= [n n])
        error('The dimension of the source power matrix does not match the source number.');
    end
    P = p;
end
end

