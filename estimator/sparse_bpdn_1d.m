function [sp, noise_var] = sparse_bpdn_1d(R, n, design, wavelength, grid_size, lambda, varargin)
%SPARSE_BPDN_1D 1D DOA estimation using basis pursuit with denosing/sparse
%recovery.
%This function utilizes 'quadprog' to solve the LASSO-like problem, and
%'SDPT3' to solve the second-order cone problem.
%Syntax:
%   sp = SPARSE_BPDN_1D(R, n, design, wavelength, grid_size, lambda, ...);
%   sp = SPARSE_BPDN_1D(R, n, f_steering, [], grid_size, lambda, ...);
%Inputs:
%   R - Sample covariance matrix.
%   design - Array design. Can also be a function handle that generates
%            a steering matrix. This function must take two arguments,
%            wavelength and the doa vector.
%   n - Number of sources.
%   wavelength - Wavelength. If design is set to a function handle, this
%                parameter must be set to [].
%   grid_size - Number of grid points used.
%   lambda - When 'Formulation' is set to 'Penalized', this is the
%            regularization parameter of the l1-norm term. When
%            'Formulation' is set to 'ConstrainedL1' or 'ConstrainedL2',
%            this is the upper bound of the l1 or l2 error.
%   ... - Options:
%           'Unit' - Can be 'radian', 'degree', or 'sin'. Default value is
%                   'radian'.
%           'Formulation' - 'Penalized' or 'ConstrainedL1' or
%               'ConstrainedL2'. The 'Penalized' formulation is give by
%                   min 0.5\|Ax - b\|_2^2 + \lambda\|x\|_1
%                   s.t. x \succeq 0
%               The 'ConstrainedL1' formulation is given by
%                   min \|Ax - b\|_2^2
%                   s.t. \|x\|_1 \leq \lambda
%                        x \succeq 0
%               while the 'ConstrainedL2' formulation is given by
%                   min \|x\|_1
%                   s.t. \|Ax - b\|_2 \leq \lambda
%                        x \succeq 0
%               Default setting is 'Penalized'.
%           'NoiseVariance' - If present, the noise variance is assumed
%               known and set to the provided value, and the objective
%               function will utilize this fact. This option is invalid
%               when the 'Formulation' is set to 'ConstrainedL2'.
%           'Verbose' - If set to true, will display detailed solver
%                       outputs.
%Output:
%   sp - Spectrum structure with the following fields:
%           x - An 1 x grid_size vector.
%           y - An 1 x grid_size vector. Calling `plot(x, y)` will plot the
%               spectrum.
%           x_est - An 1 x n vector storing the estimated DOAs.
%           x_unit - The same as the unit specified by 'Unit'.
%           resolved - True if the number of peaks in the spectrum is
%                      greater or equal to the number of sources.
%           discrete - Constant value true.
%   noise_var - If noise variance is not specified, and the formulation is
%               not set to 'ConstrainedL2', returns the estimated
%               noise variance.
unit = 'radian';
is_noise_var_known = false;
noise_var = -1;
use_constrained_formulation = false;
upper_bound_l2_norm = false;
verbose = false;
for ii = 1:2:nargin-6
    option_name = varargin{ii};
    option_value = varargin{ii+1};
    switch lower(option_name)
        case 'unit'
            unit = option_value;
        case 'formulation'
            switch lower(option_value)
                case 'penalized'
                    use_constrained_formulation = false;
                    upper_bound_l2_norm = false;
                case 'constrainedl1'
                    use_constrained_formulation = true;
                    upper_bound_l2_norm = false;
                case 'constrainedl2'
                    use_constrained_formulation = true;
                    upper_bound_l2_norm = true;
                otherwise
                    error('Unknown formulation ''%s''.', option_value);
            end
        case 'noisevariance'
            is_noise_var_known = true;
            noise_var = option_value;
        case 'verbose'
            verbose = option_value;
        otherwise
            error('Unknown option ''%s''.', option_name);
    end
end
if upper_bound_l2_norm && is_noise_var_known
    warning('Specified noise variance will be ignore when the problem formulation is set to ''ConstrainedL2''.');
end
% discretize and create the corresponding steering matrix
[doa_grid, doa_grid_display, ~] = default_doa_grid(grid_size, unit, 1);
if ishandle(design)
    A = design(wavelength, doa_grid);
else
    A = steering_matrix(design, wavelength, doa_grid);
end
% preparing the objective function
[m, ~] = size(A);
if is_noise_var_known && ~upper_bound_l2_norm
    % remove noise variance from the main diagonal if known
    R = R - noise_var*eye(m);
end
r = R(:);
Phi = khatri_rao(conj(A), A);
if ~is_noise_var_known && ~upper_bound_l2_norm
    % is noise variance is unknown, include it as a unknown parameter
    % we need to add a column here
    i = reshape(eye(m), [], 1);
    Phi = [Phi i];
end
Phi_full = [real(Phi); imag(Phi)];
r_full = [real(r); imag(r)];
[dim_s, dim_x] = size(Phi_full);
% solve
if use_constrained_formulation && upper_bound_l2_norm
    % use SDPT3
    % covert it to a SOCP problem
    check_opt_solver('sdpt3');
    blk{1,1} = 'q';
    blk{1,2} = dim_s + 1;
    At{1} = -speye(dim_s + 1, dim_s + 1);
    At{1}(1,1) = 1;
    C{1} = sparse(dim_s + 1, 1);
    blk{2,1} = 'l';
    blk{2,2} = dim_x;
    At{2} = [zeros(1, dim_x); Phi_full]';
    C{2} = ones(dim_x, 1);
    sqlp_options = sqlparameters;
    if ~verbose
        sqlp_options.printlevel = 0;
    end
    [~,x,~,~,info,~] = sqlp(blk, At, C, [lambda; r_full], sqlp_options);
    if info.termcode ~= 0
        warning('Infeasibility detected by SDPT3. Result may not be correct.');
    end
    x = x{2};
else
    % use quadprog
    check_opt_solver('quadprog');
    if ~verbose
        if exist('mskoptimset.m', 'file')
            qp_options = mskoptimset('Diagnostics', 'off');
        else
            qp_options = optimoptions(@quadprog, 'Display', 'none');
        end
    else
        qp_options = [];
    end
    H = Phi_full'*Phi_full;
    H = 0.5*(H + H');
    % force PSD
    H = H + eye(size(H)) * (H(1,1) / 1e9); % this should have minimal impact
    if use_constrained_formulation
        % constrained l1
        x = quadprog(...
            H, ...
            -r_full'*Phi_full, ...
            ones(1, dim_x), lambda, ...
            [], [], ...,
            zeros(dim_x, 1), [], ...
            [], qp_options);
    else
        % penalized l1
        x = quadprog(...
            H, ...
            lambda*ones(1, dim_x) - r_full'*Phi_full, ...
            [], [], ...
            [], [], ...
            zeros(dim_x, 1), [], ...
            [], qp_options);
    end
end
if isempty(x)
    x = zeros(dim_x, 1);
end
if ~is_noise_var_known
    if upper_bound_l2_norm
        % noise variance is not specified for the 'ConstrainedL2' case
        noise_var = lambda / sqrt(m);
    else
        noise_var = x(end);
        x = x(1:end-1);
    end
end
% find n DOAs
[x_est_intl, ~, resolved] = find_doa_from_spectrum_1d(doa_grid_display, x, n);
% return
sp = struct();
sp.x = doa_grid_display;
sp.x_est = x_est_intl;
sp.x_unit = unit;
sp.y = x';
sp.resolved = resolved;
sp.discrete = true;
end

