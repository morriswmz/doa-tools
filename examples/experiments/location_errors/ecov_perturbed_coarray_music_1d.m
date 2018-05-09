function C = ecov_perturbed_coarray_music_1d(design, wavelength, doas, p, noise_var, snapshot_count, error_stat, mode)
%ECOV_PERTURBED_COARRAY_MUSIC_1D Asymptotic covariance matrix of the
%estimation errors of the coarray-based MUSIC algorithm (SS-MUSIC and
%DA-MUSIC) in the presence of small model errors.
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
%   error_stat - A struct describing the statistical properties of the
%       model errors. The following fields are supported.
%           'PositionErrorMask' - (Optional) A boolean matrix or vector
%               specifying which position error parameter is present.
%               (1) If an 1xM vector is provided, only perturbations along
%               the x-axis (along the array) will be considered. Set the
%               m-th element of this vector to true if the m-th element is
%               perturbed.
%               (2) If a 2xM vector is provided, perturbations along both
%               the x-axis and the y-axis (perpendicular to the array) will
%               be considered. The first row corresponds to the x-axis and
%               the second row corresponds to the y-axis. For example,
%               setting the (2,m)-th element of the matrix to true if the
%               m-th element is perturbed along the y-axis.
%               (3) If not specified, the following default value will be
%               used (assuming the first element as the reference element):
%                   [false true true ... true;...
%                    false true true ... true]
%            'PositionErrorCov' - The covariance matrix of the position
%               errors. The dimension of this matrix must match the number
%               of true values in PositionErrorMask.
%   mode - 'Full' or 'DiagonalsOnly'. Default value if 'Full'.
%Output:
%   C - If mode is set to 'Full', returns the full asymptotic error
%       covariance matrix. If mode is set to 'DiagonalsOnly', will only
%       return a vector of the diagonals. This is useful when you only need
%       the asymptotic variances for each DOA.
if design.dim ~= 1
    error('1D array expected.');
end
if nargin <= 6
    mode = 'full';
    if nargin <= 5
        snapshot_count = 1;
    end
end
% Fall back to the perturbation free version if no model errors are present.
if isempty(error_stat) || isempty(fieldnames(error_stat))
    C = ecov_coarray_music_1d(design, wavelength, doas, p, noise_var, snapshot_count, mode);
    return;
end

m = design.element_count;
k = length(doas);
p = unify_source_power_vector(p, k);
% Generate the selection matrix.
F = coarray_selection_matrix_1d(design);
% Check source number.
m_v = (size(F,1)+1)/2;
if k >= m_v
    error('Too many sources.');
end

% Compute analytical MSE
% ======================

% Nominal R
design0 = design;
design0.mutual_coupling = [];
design0.position_errors = [];
design0.gain_errors = [];
design0.phase_errors = [];
A0 = steering_matrix(design0, wavelength, doas);
R0 = A0 * bsxfun(@times, p, A0') + eye(m) * noise_var;

% Coarray
arr_virtual = ula_1d(m_v, wavelength/2);
[Av, DAv] = steering_matrix(arr_virtual, wavelength, doas);
Rss = Av * bsxfun(@times, p, Av') + eye(m_v) * noise_var;
G = zeros(m_v*m_v, 2*m_v-1);

% Concatenated subarray selection matrix
for ii = 1:m_v
    G(((ii-1)*m_v+1):ii*m_v, :) = [zeros(m_v, m_v-ii) eye(m_v) zeros(m_v, ii-1)];
end
% Eigendecomposition of the ideal augmented covariance matrix
[E, ~] = eig(0.5*(Rss + Rss'), 'vector');
En = E(:,1:end-k);
% Evalute xi_k / p_k
Xi_g = zeros(m*m, k);
gammas = zeros(k, 1);
pinv_Av = pinv(Av);
for kk = 1:k
    d_ak = DAv(:,kk);
    alpha_k = -pinv_Av(kk,:).';
    beta_k = En*En'*d_ak;
    gammas(kk) = real((En'*d_ak)'*(En'*d_ak));
    Xi_g(:,kk) = F'*G'*kron(beta_k, alpha_k) / gammas(kk) / p(kk);
end
% Add model errors' contributions
% Note: Only sensor location errors are implemented.
n_perturb = 0; % number of types of perturbations
B = cell(3,1); % linearized perturbation effect matrix, \Delta r = B\delta
S = cell(3,1); % covariance matrix of the perturbation parameters
if isfield(error_stat, 'PositionErrorCov')
    if isfield(error_stat, 'PositionErrorMask')
        [nr_mask, nc_mask] = size(error_stat.PositionErrorMask);
        % dimension check
        if nc_mask ~= m
            error('The number of columns of ''PositionErrorMask'' does not match the number of array elements.');
        end
        if nr_mask < 1 || nr_mask > 2
            error('The number of rows of ''PositionErrorMask'' should be either 1 or 2.');
        end
        % convert
        if nr_mask == 1
            cur_mask = [error_stat.PositionErrorMask false(1, m)];
        else
            cur_mask = [error_stat.PositionErrorMask(1,:) error_stat.PositionErrorMask(2,:)];
        end
    else
        % assuming the first element as the reference element
        cur_mask = true(1,2*m);
        cur_mask([1 m+1]) = false;
    end
    n_perturb = n_perturb + 1;
    n_pos_err = length(find(cur_mask));
    cur_S = error_stat.PositionErrorCov;
    if isscalar(cur_S)
        % scalar case
        cur_S = cur_c * eye(n_pos_err);
    else
        if any(n_pos_err ~= size(error_stat.PositionErrorCov))
            error('The dimension %d of ''PositionErrorCov'' does not match the number of position error parameters %d.',...
                size(error_stat.PositionErrorCov, 1), n_pos_err);
        end
    end
    S{n_perturb} = cur_S;
    
    % alternative way of constructing B
%     Ts = A*bsxfun(@times, sin(doas') .* p, A');
%     Tc = A*bsxfun(@times, cos(doas') .* p, A');
%     Ms1 = zeros(m*m, m);
%     Ms2 = zeros(m*m, m);
%     Mc1 = zeros(m*m, m);
%     Mc2 = zeros(m*m, m);
%     for ii = 1:m
%         Ms1(m*(ii-1)+1:m*ii,ii) = Ts(:,ii);
%         Ms2(m*(ii-1)+1:m*ii,:) = diag(Ts(:,ii));
%         Mc1(m*(ii-1)+1:m*ii,ii) = Tc(:,ii);
%         Mc2(m*(ii-1)+1:m*ii,:) = diag(Tc(:,ii));
%     end
%     cur_B = [(Ms1 - Ms2) (Mc1 - Mc2)];
%     cur_B = (2j * pi / wavelength) * cur_B(:,cur_mask);
%     B{n_perturb} = cur_B;
    
    % Using equation (21a) (21b)
    Au = 1j*2*pi / wavelength * A0 * diag(sin(doas));
    Av = 1j*2*pi / wavelength * A0 * diag(cos(doas));
    P = diag(p);
    Bu = khatri_rao(eye(m), A0 * P * Au') + khatri_rao(conj(A0 * P * Au'), eye(m));
    Bv = khatri_rao(eye(m), A0 * P * Av') + khatri_rao(conj(A0 * P * Av'), eye(m));
    cur_B = [Bu Bv];
    B{n_perturb} = cur_B(:, cur_mask);
end

% Evaluate cov
switch lower(mode)
    case 'diagonalsonly'
        C = zeros(k,1);
        for kk = 1:k
            C(kk) = real(Xi_g(:,kk)' * kron(R0, R0.') * Xi_g(:,kk)) / snapshot_count;
            for pp = 1:n_perturb
                T = real(Xi_g(:,kk).' * B{pp});
                C(kk) = C(kk) + T * S{pp} * T';
            end
        end
    case 'full'
        C = real(Xi_g' * kron(R0, R0.') * Xi_g / snapshot_count);
        for pp = 1:n_perturb
            T = real(Xi_g.' * B{pp});
            C = C + T * S{pp} * T';
        end
    otherwise
        error('Invalid mode.');
end
end

