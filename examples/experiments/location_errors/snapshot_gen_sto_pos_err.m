function [X, R, S] = snapshot_gen_sto_pos_err(design, doas, wavelength, t, ncov, scov, pos_err_info)
%SNAPSHOT_GEN_STO_POS_ERR Generates snapshots for the stochastic model
%with 2D sensor location errors changes every snapshot.
%Syntax:
%   X = SNAPSHOT_GEN_STO_POS_ERR_2D(design, doas, wavelength[, t, ncov, scov, pos_err_info]);
%   [X, R] = SNAPSHOT_GEN_STO_POS_ERR_2D(design, doas, wavelength[, t, ncov, scov, pos_err_info]);
%   [X, R, S] = SNAPSHOT_GEN_STO_POS_ERR_2D(design, doas, wavelength[, t, ncov, scov, pos_err_info]);
%Inputs:
%   design - Array design.
%   doas - DOA vector. For 2D DOAs, each column represents a DOA pair.
%   wavelength - Wavelength.
%   t - Number of snapshots.
%   ncov - Covariance matrix of the additive complex circular-symmetric
%          Gaussian noise. Can be a scalar, vector (for uncorrelated noise
%          with different powers), or a matrix.
%   scov - Covariance matrix of the source signals. Can be a scalar, vector
%          (for uncorrelated sources with different powers), or a matrix.
%   pos_err_info - A struct containing the information about position
%                  errors:
%                   'std' - Standard deviation.
%                   'mask' - Specifies which sensor has location errors.
%                   'reference' - 'Absolute' or 'First'.
%Outputs:
%   X - Snapshots, where each columns is a single snapshot.
%   R - Sample covariance matrix (averaged by the number of snapshots).
%   S - A source_count x snapshot_count matrix consists of source signal
%       vectors.
if nargin <= 6
    pos_err_info = struct('std', 0, 'mask', true(design.element_count, 1), ...
        'reference', 'Absolute');
end
if nargin <= 5
    scov = 1;
end
if nargin <= 4
    ncov = 1;
end
if nargin <= 3
    t = 1;
end
m = design.element_count;
k = length(doas);
S_internal = gen_ccsg(k, t, scov);
N = gen_ccsg(m, t, ncov);
X = zeros(m, t);
ind = find(pos_err_info.mask);
n_err = length(ind);
% lazy comparison here
ref_first = strcmpi(pos_err_info.reference, 'first');
for tt = 1:t
    % generate position errors
    pos_err = zeros(2, m);
    pos_err(:,ind) = pos_err_info.std * randn(2, n_err);
    if ref_first
        pos_err = pos_err - pos_err(:,1);
    end
    design.position_errors = pos_err;
    A = steering_matrix(design, wavelength, doas);
    X(:,tt) = A * S_internal(:,tt) + N(:,tt);
end
if nargout >= 2
    R = (X*X')/t;
    if nargout == 3
        S = S_internal;
    end
end
end

function X = gen_ccsg(m, n, cov)
X0 = randn(m, n) + 1j * randn(m, n);
if isscalar(cov)
    X = sqrt(cov/2) * X0;
elseif isvector(cov)
    X = bsxfun(@times, X0, cov(:));
else
    C = sqrtm(cov/2);
    X = C*X0;
end
end