function sp = music_1d(R, n, design, wavelength, grid_size, varargin)
%MUSIC_1D 1D MUSIC.
%Syntax:
%   sp = music_1d(R, n, design, wavelength, grid_size, ...);
%   sp = music_1d(R, n, f_steering, [], grid_size, ...);
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
m = size(R, 1);
if n >= m
    error('Too many sources.');
end
need_deg_conv = false;
switch lower(unit)
    case 'radian'
        doa_grid = linspace(-pi/2, pi/2, grid_size);
    case 'degree'
        doa_grid = linspace(-pi/2, pi/2, grid_size);
        need_deg_conv = true;
    case 'sin'
        doa_grid = asin(linspace(-1, 1, grid_size));
    otherwise
        error('Unknown unit "%s".', unit);
end
if ishandle(design) && isempty(wavelength)
    A = design(wavelength, doa_grid);
else
    A = steering_matrix(design, wavelength, doa_grid);
end
% find noise subspace
[U, D] = eig(0.5*(R + R'));
% possible asymmetry due to floating point error
if ~isreal(D)
    eig_values = abs(diag(D));
    [~, I] = sort(eig_values);
    Un = U(:, I(1:end-n));
else
    Un = U(:, 1:end-n);
end
% compute spectrum
sp_intl = Un'*A;
sp_intl = sum(real(sp_intl).^2 + imag(sp_intl).^2, 1)';
sp_intl = 1./sp_intl;
[peaks, idx_est] = findpeaks(sp_intl);
if length(idx_est) == n
    % juse fine
    x_est_intl = doa_grid(idx_est);
    resolved = true;
elseif length(idx_est) > n
    % more than expected, choose largest ones
    [~, sort_idx] = sort(peaks);
    x_est_intl = doa_grid(idx_est(sort_idx(end-n+1:end)));
    x_est_intl = sort(x_est_intl);
    resolved = true;
else
    % failed to detect
    x_est_intl = [];
    resolved = false;
end
% return
sp = struct();
sp.x_est = x_est_intl';
if need_deg_conv
    sp.x = rad2deg(doa_grid)';
else
    sp.x = doa_grid';
end
sp.x_unit = unit;
sp.y = sp_intl;
sp.resolved = resolved;
sp.discrete = false;
end

