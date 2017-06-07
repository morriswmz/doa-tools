function [X, R, S] = snapshot_gen_sto(design, doas, wavelength, t, ncov, scov)
%SNAPSHOT_GEN_STO Generates snapshots for the stochastic model.
%Syntax:
%   X = STO_SNAPSHOT_GEN(design, doas, wavelength[, t, ncov, scov]);
%   [X, R] = STO_SNAPSHOT_GEN(design, doas, wavelength[, t, ncov, scov]);
%   [X, R, S] = STO_SNAPSHOT_GEN(design, doas, wavelength[, t, ncov, scov]);
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
%Outputs:
%   X - Snapshots, where each columns is a single snapshot.
%   R - Sample covariance matrix (averaged by the number of snapshots).
%   S - A source_count x snapshot_count matrix consists of source signal
%       vectors.
if nargin <= 5
    scov = 1;
end
if nargin <= 4
    ncov = 1;
end
if nargin <= 3
    t = 1;
end
A = steering_matrix(design, wavelength, doas);
[m, k] = size(A);
S_internal = gen_ccsg(k, t, scov);
X = A * S_internal + gen_ccsg(m, t, ncov);
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