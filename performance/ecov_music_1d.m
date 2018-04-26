function C = ecov_music_1d(design, wavelength, doas, P, noise_var, snapshot_count)
%ECOV_MUSIC_1D Asymptotic covariance matrix of the estimation errors
%of the classical MUSIC algorithm.
%Syntax:
%   C = ECOV_MUSIC_1D(design, wavelength, doas, P, noise_var[, snapshot_count])
%Inputs:
%   design - Array design.
%   wavelength - Wavelength.
%   doas - DOA vector in radians.
%   P - Source covariance matrix. If all sources are uncorrelated and
%       shares the same power, you can just pass in a scalar. If all
%       sources are uncorrelated but have different powers, you can just
%       pass in a vector.
%   noise_var - Noise power.
%   snapshot_count - (Optional) number of snapshots. Default is one.
%Output:
%   C - Asymptotic error covariance matrix.
%Reference:
%   [1] P. Stoica and A. Nehorai, "MUSIC, maximum likelihood, and
%       Cramer-Rao bound: further results and comparisons," IEEE
%       Transactions on Acoustics, Speech and Signal Processing, vol. 38,
%       no. 12, pp. 2140-2150, Dec. 1990.
%   [2] P. Stoica and A. Nehorai, "MUSIC, maximum likelihood, and
%       Cramer-Rao bound," IEEE Transactions on Acoustics, Speech and
%       Signal Processing, vol. 37, no. 5, pp. 720-741, May 1989.
if design.dim ~= 1
    error('1D array expected.');
end
if nargin <= 5
    snapshot_count = 1;
end
m = design.element_count;
k = length(doas);
P = unify_source_power_matrix(P, k);
[A, D] = steering_matrix(design, wavelength, doas);
H = D'*(eye(m) - A/(A'*A)*A')*D;
B = P\(P + noise_var*eye(k)/(A'*A))/P;
h = diag(H);
C = (noise_var/2/snapshot_count) * (real(H .* B) .* (h*h'));
end