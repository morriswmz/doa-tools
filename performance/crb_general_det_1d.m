function CRB = crb_general_det_1d(design, wavelength, doas, P_est, noise_var, snapshot_count)
%CRB_GENERAL_DET_1D CRB for general 1D arrays based on the deterministic
%(conditional) model, in radians.
%Inputs:
%   design - Array design.
%   wavelength - Wavelength.
%   doas - DOA vector in radians.
%   P_est - Estimated source covariance matrix from deterministic source
%           signals.
%   noise_var - Noise power.
%   snapshot_count - (Optional) number of snapshots. Default is one.
%Reference:
%   [1] P. Stoica and A. Nehorai, "Performance study of conditional and
%       unconditional direction-of-arrival estimation," IEEE Transactions
%       on Acoustics, Speech and Signal Processing, vol. 38, no. 10,
%       pp. 1783–1795, Oct. 1990.
if design.dim ~= 1
    error('1D array expected.');
end
if nargin <= 5
    snapshot_count = 1;
end
[A, D] = steering_matrix(design, wavelength, doas);
[m, k] = size(A);
H = D'*(eye(m) - A/(A'*A)*A')*D;
CRB = real(H .* P_est.');
CRB = eye(k) / CRB * (noise_var / snapshot_count / 2);
end

