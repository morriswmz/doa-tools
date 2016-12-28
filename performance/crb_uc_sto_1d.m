function [CRB, FIM] = crb_uc_sto_1d(design, wavelength, doas, p, noise_var, snapshot_count)
%CRB_UC_STO_1D CRB Stochastic CRB for general 1D arrays with the assumption
%that the sources are uncorrelated.
if design.dim ~= 1
    error('1D array expected.');
end
if nargin <= 5
    snapshot_count = 1;
end
m = design.element_count;
k = length(doas);
if ~isscalar(p) 
    if length(p) ~= k
        error('The length of the source power vector does not match number of sources.');
    end
    p = reshape(p, [], 1);
else
    p = p * ones(k, 1);
end
% we need to compute each submatrix of the FIM
[A, DA] = steering_matrix(design, wavelength, doas);
R = A * bsxfun(@times, p, A') + noise_var * eye(m);
R_inv = eye(m) / R;
R_inv = 0.5 * (R_inv + R_inv');
DRD = DA' * R_inv * DA;
DRA = DA' * R_inv * A;
ARD = A' * R_inv * DA;
ARA = A' * R_inv * A;
PP = p*p';
FIM_tt = 2*real((DRD.' .* ARA + DRA.' .* ARD) .* PP);
FIM_pp = real(ARA.' .* ARA);
R_inv2 = R_inv * R_inv;
FIM_ss = real(sum(diag(R_inv2)));
FIM_tp = 2*real(DRA.' .* (bsxfun(@times, p, ARA)));
FIM_ts = 2*real(p .* sum(DA .* (R_inv2 * A), 1)');
FIM_ps = real(sum(A .* (R_inv2 * A), 1)');
FIM = [...
    FIM_tt  FIM_tp  FIM_ts; ...
    FIM_tp' FIM_pp  FIM_ps; ...
    FIM_ts' FIM_ps' FIM_ss] * snapshot_count;
CRB = eye(2*k+1) / FIM;
CRB = 0.5 * (CRB + CRB');
CRB = CRB(1:k, 1:k);
end