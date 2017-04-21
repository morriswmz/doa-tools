function C = ecov_coarray_music_1d(design, wavelength, doas, p, noise_var, snapshot_count, mode)
%ECOV_COARRAY_MUSIC_1D Asymptotic covariance matrix of the estimation
%errors of the coarray-based MUSIC algorithm (SS-MUSIC and DA-MUSIC).
%Syntax:
%   C = ECOV_COARRAY_MUSIC_1D(design, wavelength, doas, P, noise_var, snapshot_count);
%   C = ECOV_COARRAY_MUSIC_1D(design, wavelength, doas, P, noise_var, snapshot_count, 'DiagonalsOnly');
%Inputs:
%   design - Array design.
%   wavelength - Wavelength.
%   doas - DOA vector in radians.
%   p - Source power vector (diagonals of the source covariance matrix).
%       If all sources share the same power, you can just pass in a scalar.
%       You can also pass in a diagonal matrix, whose diagonals will be
%       extracted automatically.
%   noise_var - Noise power.
%   snapshot_count - (Optional) number of snapshots. Default is one.
%Output:
%   C - If mode is set to 'Full', returns the full asymptotic error
%       covariance matrix. If mode is set to 'DiagonalsOnly', will only
%       return a vector of the diagonals. This is useful when you only need
%       the asymptotic variances.
if design.dim ~= 1
    error('1D array expected.');
end
if nargin <= 6
    mode = 'full';
    if nargin <= 5
        snapshot_count = 1;
    end
end
m = design.element_count;
k = length(doas);
p = unify_source_power_vector(p, k);
% generate selection matrix
F = coarray_selection_matrix_1d(design);
% check source number
m_v = (size(F,1)+1)/2;
if k >= m_v
    error('Too many sources.');
end
% compute analytical MSE here
A = steering_matrix(design, wavelength, doas);
R_ideal = A * bsxfun(@times, p, A') + eye(m) * noise_var;
arr_virtual = ula_1d(m_v, wavelength/2);
[Av, DAv] = steering_matrix(arr_virtual, wavelength, doas);
Rss_ideal = Av * bsxfun(@times, p, Av') + eye(m_v) * noise_var;
G = zeros(m_v*m_v, 2*m_v-1);
% concatenated subarray selection matrix
for ii = 1:m_v
    G(((ii-1)*m_v+1):ii*m_v, :) = [zeros(m_v, m_v-ii) eye(m_v) zeros(m_v, ii-1)];
end
% eigendecomposition of the ideal augmented covariance matrix
[E, ~] = eig(0.5*(Rss_ideal + Rss_ideal'), 'vector');
En = E(:,1:end-k);
% evalute xi_k
Xi_g = zeros(m*m, k);
gammas = zeros(k, 1);
pinv_Av = pinv(Av);
for kk = 1:k
    d_ak = DAv(:,kk);
    alpha_k = -pinv_Av(kk,:).' / p(kk);
    beta_k = En*En'*d_ak;
    gammas(kk) = real((En'*d_ak)'*(En'*d_ak));
    Xi_g(:,kk) = F'*G'*kron(beta_k, alpha_k) / gammas(kk);
end
% evaluate cov
switch lower(mode)
    case 'diagonalsonly'
        C = zeros(k,1);
        for kk = 1:k
            C(kk) = real(Xi_g(:,kk)'*kron(R_ideal, R_ideal.')*Xi_g(:,kk))/snapshot_count;
        end
    case 'full'
        C = real(Xi_g'*kron(R_ideal, R_ideal.')*Xi_g/snapshot_count);
    otherwise
        error('Invalid mode.');
end
end

