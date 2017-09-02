function [X, R, S] = snapshot_gen_sym(design, doas, wavelength, stype, t, ncov, sp)
%SNAPSHOT_GEN_SYM Generates snapshots from symbols.
%Syntax:
%   X = SNAPSHOT_GEN_SYM(design, doas, wavelength, stype[, t, ncov, sp]);
%   [X, R] = SNAPSHOT_GEN_SYM(design, doas, wavelength, stype[, t, ncov, sp]);
%   [X, R, S] = SNAPSHOT_GEN_SYM(design, doas, wavelength, stype[, t, ncov, sp]);
%Inputs:
%   design - Array design.
%   doas - DOA vector. For 2D DOAs, each column represents a DOA pair.
%   wavelength - Wavelength.
%   stype - Symbol type.
%   t - Number of snapshots.
%   ncov - Covariance matrix of the additive complex circular-symmetric
%          Gaussian noise. Can be a scalar, vector (for uncorrelated noise
%          with different powers), or a matrix.
%   sp - Sources powers. 
%Outputs:
%   X - Snapshots, where each columns is a single snapshot.
%   R - Sample covariance matrix (averaged by the number of snapshots).
%   S - A source_count x snapshot_count matrix consists of source signal
%       vectors.
if nargin <= 6
    sp = 1;
end
if nargin <= 5
    ncov = 1;
end
if nargin <= 4
    t = 1;
end
A = steering_matrix(design, wavelength, doas);
[m, k] = size(A);
S_internal = gen_symbols(k, t, stype, sp);
X = A * S_internal + gen_ccsg(m, t, ncov);
if nargout >= 2
    R = (X*X')/t;
    if nargout == 3
        S = S_internal;
    end
end
end

function S = gen_symbols(m, n, stype, sp)
switch lower(stype)
    case 'bpsk'
        symbols = [-1 1];
    case 'qpsk'
        symbols = exp(1j*[1 3 5 7]/4*pi);
    otherwise
        error('Unsupported symbol type.')
end
S = symbols(randi([1 length(symbols)], m, n));
if isscalar(sp)
    S = sqrt(sp) * S;
else
    S = bsxfun(@times, sqrt(sp(:)), S);
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