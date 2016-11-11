function CRB = crb_general_sto_1d(design, wavelength, doas, P, noise_var, snapshot_count)
%CRB_GENERAL_STO_1D CRB for general 1D arrays based on the stochastic
%(unconditional) model, in radians.
%Inputs:
%   design - Array design.
%   wavelength - Wavelength.
%   doas - DOA vector.
%   P - Source covariance matrix.
%   noise_var - Noise power.
%   snapshot_count - (Optional) number of snapshots. Default is one.
%Reference:
% [1] P. Stoica and A. Nehorai, "Performance study of conditional and
%     unconditional direction-of-arrival estimation," IEEE Transactions on
%     Acoustics, Speech and Signal Processing, vol. 38, no. 10,
%     pp. 1783–1795, Oct. 1990.
if design.dim ~= 1
    error('1D array expected.');
end
if nargin <= 5
    snapshot_count = 1;
end
A = steering_matrix(design, wavelength, doas);
[m, k] = size(A);
R = A*P*A' + noise_var * eye(m);
D = zeros(m, k);
% compute derivative w.r.t. DOAs
for kk = 1:k
    D(:,kk) = (2j*pi/wavelength*design.element_positions(:)*cos(doas(kk))) ...
                .* A(:,kk);
end
H = D'*(eye(m) - A/(A'*A)*A')*D;
CRB = real(H .* ((P*A'/R*A*P).'));
CRB = eye(k) / CRB * (noise_var / snapshot_count / 2);
end

