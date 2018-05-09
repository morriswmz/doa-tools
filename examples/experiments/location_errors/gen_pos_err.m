function [err, C] = gen_pos_err(std, n, type, absolute)
%GEN_POS_ERR Generates position errors. We first generate i.i.d. position errors
%for each sensor. We use the first sensor as the reference sensor and substract
%the the position error of the first sensor from the remaining sensors.
%Therefore, the resulting covariance matrix is not diagonal.
%Inputs:
%   std - Standard deviation of the position errors.
%   n - Number of sensors.
%   type - 'Gaussian' or 'Uniform'.
%   absolute - If set to true, will not use the first sensor as the reference
%              sensor and the position errors will be i.i.d.
%Outputs:
%   err - 2xn vector of position errors.
%   C - Covariance matrix. (2n-2) x (2n-2) if absolute is false, (2n) x (2n) if
%       absolute is true.
if nargin < 4
    absolute = false;
end
switch lower(type)
    case 'gaussian'
        err = randn(2, n)*std;
        if nargout > 1
            if absolute
                C = std^2*eye(n + n);
            else
                C = std^2*(eye(n - 1) + ones(n - 1));
                C = blkdiag(C, C);
            end
        end
    case 'uniform'
        b = sqrt(3)*std;
        err = unifrnd(-b, b, 2, n);
        if nargout > 1
            if absolute
                C = std^2*eye(n + n);
            else
                C = std^2*(eye(n - 1) + ones(n - 1));
                C = blkdiag(C, C);
            end
        end
    otherwise
        error('Unknow type "%s".', type);
end
if ~absolute
    err(1,2:end) = err(1,2:end) - err(1,1);
    err(2,2:end) = err(2,2:end) - err(2,1);
    err(1,1) = 0;
    err(2,1) = 0;
end
end

