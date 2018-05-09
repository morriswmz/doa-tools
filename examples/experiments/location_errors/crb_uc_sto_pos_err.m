function [CRB, FIM] = crb_uc_sto_pos_err(design, wavelength, doas, p, noise_var, snapshot_count, pos_err_mask)
%CRB_UC_STO_POS_ERR Stochastic CRB for joint estimation of DOA parameters and
%the deterministic location errors with the assumption that the sources are
%uncorrelated.
%Inputs:
%   design - Array design.
%   wavelength - Wavelength.
%   doas - DOA vector.
%   p - Source power vector (diagonals of the source covariance matrix).
%       If all sources share the same power, you can just pass in a scalar.
%       You can also pass in a diagonal matrix, whose diagonals will be
%       extracted automatically.
%   noise_var - Noise power.
%   snapshot_count - (Optional) number of snapshots. Default is one.
%   pos_err_mask - (Optional) A 2xn or 1xn boolean mask specifying the position
%                  error parameters to be considered. If omitted, the masks will
%                  be constructed from the non-zero elements from
%                  design.position_errors.
%Outputs:
%   CRB - The CRB matrix.
%   FIM - The full FIM.
if design.dim ~= 1
    error('1D array expected.');
end
if nargin <= 5
    snapshot_count = 1;
end
m = design.element_count;
% Parse position error masks.
if nargin <= 6
    % No mask given. Generate from the position errors in the given design.
    [mask_u, mask_v] = get_pos_err_mask_from_design(design);
else
    if size(pos_err_mask, 1) == 1
        mask_u = pos_err_mask;
        mask_v = true(1, m);
    else
        mask_u = pos_err_mask(1,:);
        mask_v = pos_err_mask(2,:);
    end
end
% Fall back to the original CRB formula if no position errors are present.
if ~any(mask_u) && ~any(mask_v)
    [CRB, FIM] = crb_uc_sto_1d(design, wavelength, doas, p, noise_var, snapshot_count);
    return;
end
% Calculate the CRB according to (33)
k = length(doas);
p = unify_source_power_vector(p, k);
compute_u = any(mask_u);
compute_v = any(mask_v);
I = eye(m);
L1 = I(:, mask_u);
L2 = I(:, mask_v);
% Note that we need to use the perturbed steering matrix here.
[A, DA] = steering_matrix(design, wavelength, doas);
R = A * bsxfun(@times, p, A') + noise_var * eye(m);
% Compute each blocks.
% Note: the Kronecker product may not be efficient for large arrays.
M_theta = (khatri_rao(conj(DA), A) + khatri_rao(conj(A), DA)) * diag(p);
M_p = khatri_rao(conj(A), A);
if compute_u
    Au = 1j*2*pi/wavelength * A * diag(sin(doas));
    APAu = A * diag(p) * Au';
    M_u = khatri_rao(conj(APAu) * L1, L1) + khatri_rao(L1, APAu * L1); 
else
    M_u = zeros(m*m, 0);
end
if compute_v
    Av = 1j*2*pi/wavelength * A * diag(cos(doas));
    APAv = A * diag(p) * Av';
    M_v = khatri_rao(conj(APAv) * L2, L2) + khatri_rao(L2, APAv * L2); 
else
    M_v = zeros(m*m, 0);
end
M_n = I(:);
M = [M_theta M_p M_u M_v M_n];
FIM = M' / kron(R.', R) * M;
FIM = 0.5 * real(FIM + FIM') * snapshot_count;
CRB = eye(size(FIM)) / FIM;
CRB = CRB(1:k, 1:k);
end

