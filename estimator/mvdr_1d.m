function sp = mvdr_1d(R, n, design, wavelength, grid_size, varargin)
%MVDR_1D 1D MVDR beamforming.
%Syntax:
%   sp = MVDR_1D(R, n, design, wavelength, grid_size, ...);
%   sp = MVDR_1D(R, n, f_steering, [], grid_size, ...);
%Inputs:
%   R - Sample covariance matrix.
%   n - Number of sources.
%   design - Array design. Can also be a function handle that generates
%            a steering matrix. This function must take two arguments,
%            wavelength and the doa vector.
%   wavelength - Wavelength. If design is set to a function handle, this
%                parameter must be set to [].
%   grid_size - Number of grid points used.
%   ... - Options:
%           'Unit' - Can be 'radian', 'degree', or 'sin'. Default value is
%                   'radian'.
%Output:
%   sp - Spectrum.
unit = 'radian';
for ii = 1:2:nargin-5
    option_name = varargin{ii};
    option_value = varargin{ii+1};
    switch lower(option_name)
        case 'unit'
            unit = option_value;
        otherwise
            error('Unknown option "%s".', option_name);
    end
end
% discretize and create the corresponding steering matrix
[A, ~, doa_grid_display, ~] = default_steering_matrix_grid(design, wavelength, grid_size, unit, 1);
% compute spectrum
R_inv = eye(size(R)) / R;
sp_intl = zeros(1, grid_size);
% because we cannot determine if we can compute the square root of R^{-1}
% we will use the for-loop implementation here.
for ii = 1:grid_size
    sp_intl(ii) = real(A(:,ii)' * R_inv * A(:,ii));
end
sp_intl = 1./sp_intl;
[x_est_intl, resolved] = find_doa_est_1d(doa_grid_display, sp_intl, n);
% return
sp = struct();
sp.x = doa_grid_display;
sp.x_est = x_est_intl;
sp.x_unit = unit;
sp.y = sp_intl;
sp.resolved = resolved;
sp.discrete = false;
end


