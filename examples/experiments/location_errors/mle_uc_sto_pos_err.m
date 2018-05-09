function [doas, p, u, v, noise_var] = mle_uc_sto_pos_err(R, n_doas, design, wavelength, x0, pos_err_mask)
%MLE_UC_STO_POS_ERR Stochastic maximum-likelihood estimator with the
%uncorrelated source assumption for jointly estimate the DOA parameters and
%the sensor location errors.
%Inputs:
%   R - Sample covariance matrix.
%   n_doas - Number of sources.
%   design - Array design.
%   wavelength - Wavelength.
%   x_0 - (Optional) an vector storing the initial guess in the following order:
%         [doas; p; u; v; noise_var]. If absent, an initial guess will be
%         obtained using SS-MUSIC. The input DOAs must be in radians.
%   pos_err_mask - (Optional) A 2xn or 1xn boolean mask specifying the position
%                  error parameters to be considered. If omitted, the masks will
%                  be constructed from the non-zero elements from
%                  design.position_errors.
%Outputs:
%   doas - DOA estimates in radians.
%   p - Estimates of the source powers.
%   u - Estimates of the position errors along the x-axis.
%   v - Estimates of the position errors along the y-axis.
%   noise_var - Estimate of the noise variance.
if nargin < 5
    x0 = [];
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
% Find number of unknown location errors.
n_u = length(find(mask_u));
n_v = length(find(mask_v));
n_unknown = n_doas + n_doas + n_u + n_v + 1;
if ~isempty(x0) && length(x0) ~= n_unknown
    error('The vector of the initial guess has incorrect length.');
end
design_nominal = design;
design_nominal.position_errors = [];
% Compute the initial guess.
if isempty(x0)
    x0 = zeros(n_unknown, 1);
    % Use SS-MUSIC to obtain the initial guess of the DOAs.
    [Rss, ~] = virtual_ula_cov_1d(design, R, 'SS');
    sp = rmusic_1d(Rss, n_doas, 2 * pi * design.element_spacing / wavelength);
    x0(1:n_doas) = sp.x_est;
    % Use least squares to obtain the initial guess of the source powers
    % and noise power.
    A_est = steering_matrix(design_nominal, wavelength, sp.x_est);
    pn_est = real(linsolve([khatri_rao(conj(A_est), A_est) reshape(eye(m), [], 1)], R(:)));
    x0(n_doas+1:n_doas+n_doas) = pn_est(1:end-1);
    x0(end) = pn_est(end);
    % Position errors remain zeros.
end
% NLL wrapper - decompose x into difference components and pass them to the NLL.
f = @(x) nll(R, design_nominal, wavelength, ...
    x(1:n_doas), x(n_doas+1:n_doas+n_doas), ...
    x(n_doas+n_doas+1:n_doas+n_doas+n_u), ...
    x(n_doas+n_doas+n_u+1:n_doas+n_doas+n_u+n_v), ...
    mask_u, mask_v, x(end));
% Constraints.
lb = [ones(n_doas, 1) * (-pi/2);...
      zeros(n_doas, 1);...
      ones(n_u + n_v, 1) * (-wavelength);...
      0];
ub = [ones(n_doas, 1) * pi/2;...
      inf(n_doas, 1);...
      ones(n_u + n_v, 1) * wavelength;...
      inf];
% Call fmincon. This will take a long time.
options = optimoptions('fmincon', 'MaxFunctionEvaluations', 50000, 'Display', 'none');
sol = fmincon(f, x0, [], [], [], [], lb, ub, [], options);
% Output the solution.
doas = sol(1:n_doas);
p = sol(n_doas+1:n_doas+n_doas);
u = sol(n_doas+n_doas+1:n_doas+n_doas+n_u);
v = sol(n_doas+n_doas+n_u+1:n_doas+n_doas+n_u+n_v);
noise_var = sol(end);
end

function obj = nll(R, design, wavelength, doas, p, u, v, mask_u, mask_v, noise_var)
%NLL Negative log-likelihood function to be optimized.
pos_errs = zeros(2, design.element_count);
pos_errs(1,mask_u) = u;
pos_errs(2,mask_v) = v;
design_perturbed = design;
design_perturbed.position_errors = pos_errs;
A = steering_matrix(design_perturbed, wavelength, doas);
S = A * diag(p) * A' + noise_var*eye(design.element_count);
S = 0.5 * (S + S');
[L,p] = chol(S);
if p > 0
    obj = inf;
    return;
end
logdetS = sum(log(diag(L)));
obj = logdetS + sum(real(diag(R / S)));
end
