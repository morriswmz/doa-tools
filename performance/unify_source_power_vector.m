function p = unify_source_power_vector(P, n)
%UNIFY_SOURCE_POWER_VECTOR A utility function that unifies scalar,
%vector, or matrix inputs of the source power parameter into a vector.
%Syntax:
%   p = UNIFY_SOURCE_POWER_VECTOR(P, n);
%Inputs:
%   P - A scalar, vector or matrix.
%       If P is a scalar, returns p = P*ones(n,1).
%       If P is a vector, checks its length and returns P.
%       If P is a matrix, checks its dimension and returns diag(P).
%   n - Length of the final vector.
%Outputs:
%   P - Unified source power vector.
if ~isscalar(P)
    if isvector(P)
        if length(P) ~= n
            error('The length of the source power vector does not match number of sources.');
        end
        p = reshape(P, [], 1);
    elseif isdiag(P)
        p = diag(P);
    else
        error('The source covariance matrix should be diagonal.');
    end
else
    p = P * ones(n, 1);
end
end

