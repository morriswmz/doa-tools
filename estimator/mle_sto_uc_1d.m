function [sp, p, noise_var] = mle_sto_uc_1d(R, n, design, wavelength, x0, varargin)
%MLE_STO_UC_1D Maximum-likelihood estimator for the stochastic model with
%uncorrelated source assumption. The estimates are obtained by solving the
%following optimization problem:
%   min_{\theta, p, \sigma} logdet(S) + tr(S^{-1} R)
%   s.t. S = A(\theta) diag(p) A(\theta)^H + \sigma I,
%        -pi/2 <= \theta <= pi/2,
%        p, \sigma >= 0.
%This optimization problem is solved by MATLAB's built-in function
%`fmincon`. It is adviced to provide a good initial guess to obtain good
%results.
%Syntax:
%   sp = MLE_STO_UC_1D(R, n, design, wavelength, [x_0, ...]);
%   sp = MLE_STO_UC_1D(R, n, f_steering, [], [x_0, ...]);
%   [sp, p, noise_var] = MLE_STO_UC_1D(R, n, design, wavelength, [x_0, ...]);
%   [sp, p, noise_var] = MLE_STO_UC_1D(R, n, f_steering, [], [x_0, ...]);
%Inputs:
%   R - Sample covariance matrix.
%   n - Number of sources.
%   design - Array design. Can also be a function handle that generates
%            a steering matrix. This function must take two arguments,
%            wavelength and the doa vector.
%   wavelength - Wavelength.
%   x_0 - (Optional) an (2n+1)x1 vector storing the initial guess in the
%         following order: [doas; p; noise_var]. If absent, an initial
%         guess will be obtained using the MVDR beamformer.
%   ... - Options:
%       'Unit' - Can be 'radian', 'degree', or 'sin'. Default value is
%                'radian'.
%       'Verbose' - If set to true, will display detailed solver
%                   outputs.
%Output:
%   sp - Spectrum structure with the following fields:
%           x - An 1 x grid_size vector.
%           y - An 1 x grid_size vector. Calling `plot(x, y)` will plot the
%               spectrum.
%           x_est - An 1 x n vector storing the estimated DOAs.
%           x_unit - The same as the unit specified by 'Unit'.
%           resolved - Constant value true.
%           discrete - Constant value true.
%   p - An nx1 vector of estimated source powers.
%   noise_var - Estimated noise power.
if nargin < 5
    x0 = [];
else
    if length(x0) ~= 2*n + 1
        error('The vector of the initial guess has incorrect length.');
    end
end
unit = 'radian';
verbose = false;
for ii = 1:2:nargin - 5
    option_name = varargin{ii};
    option_value = varargin{ii + 1};
    switch lower(option_name)
        case 'unit'
            unit = option_value;
        case 'verbose'
            verbose = option_value;
        otherwise
            error('Unknown option ''%s''.', option_name);
    end
end
% prepare initial guess
if isempty(x0)
    % generate initial guess
    x0 = zeros(2*n + 1, 1);
    sp_mvdr = mvdr_1d(R, n, design, wavelength, 180);
    if sp_mvdr.resolved
        x0(1:n) = sp_mvdr.x_est;
    else
        warning('Failed to obtain a good initial guess.');
        % we just assume the doas is uniformly placed
        x0(1:n) = linspace(-pi/3, pi/3, n);
    end
    if ishandle(design)
        A = design(wavelength, x0(1:n));
    else
        A = steering_matrix(design, wavelength, x0(1:n));
    end
    B = [khatri_rao(conj(A), A) reshape(eye(size(A, 1)), [], 1)];
    z = real(B\R(:));
    z(z < 0) = 0;
    x0(n + 1:end) = z;
else
    x0 = x0(:);
    switch (unit)
        case 'degree'
            x0(1:n) = deg2rad(x0(1:n));
        case 'sin'
            x0(1:n) = asin(x0(1:n));
        case 'radian'
        otherwise
            error('Unexpected unit ''%s''.', unit);
    end
end
% prepare the optimization problem
if verbose
    options = optimoptions('fmincon', 'Display', 'iter');
else
    options = optimoptions('fmincon', 'Display', 'off');
end
lb = zeros(2*n + 1, 1);
lb(1:n) = -pi/2;
ub = inf(2*n + 1, 1);
ub(1:n) = pi/2;
% solve it
x_opt = fmincon(@(x) nll(R, design, wavelength, x(1:n)', x(n + 1:n + n), x(end)), ...
    x0, [], [], [], [], lb, ub, [], options);
doas = sort(x_opt(1:n));
p = x_opt(n + 1:n + n);
noise_var = x_opt(end);
switch (unit)
    case 'degree'
        doas = rad2deg(doas);
    case 'sin'
        doas = sin(doas);
    case 'radian'
    otherwise
        error('Unexpected unit ''%s''.', unit);
end
% store results
sp = struct();
sp.x_est = doas';
sp.x = sp.x_est;
sp.x_unit = unit;
sp.y = ones(1, n);
sp.resolved = true;
sp.discrete = true;
end

function obj = nll(R, design, wavelength, doas, p, noise_var)
% negative log-likelihood function.
if ishandle(design)
    A = design(wavelength, doas);
else
    A = steering_matrix(design, wavelength, doas);
end
S = A*bsxfun(@times, p, A') + eye(size(A, 1))*noise_var;
S = 0.5*(S + S');
[L,p] = chol(S);
if p > 0
    obj = inf;
    return;
end
logdetS = sum(log(diag(L)));
obj = logdetS + sum(real(diag(R/S)));
end
