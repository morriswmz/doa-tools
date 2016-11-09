function sp = rmusic_1d(R, n, k, varargin)
%RMUSIC_1D 1D root MUSIC.
%Syntax:
%   sp = RMUSIC_1D(R, n, k);
%Inputs:
%   R - Sample covariance matrix.
%   n - Number of sources.
%   k - 2*pi*inter_element_spacing/wavelength
%   ... - Options:
%       'Unit' - Can be 'radian', 'degree', or 'sin'. Default value is
%                'radian'.
%Output:
%   sp - Spectrum struct.
unit = 'radian';
for ii = 1:2:nargin-3
    option_name = varargin{ii};
    option_value = varargin{ii+1};
    switch lower(option_name)
        case 'unit'
            unit = option_value;
        otherwise
            error('Unknow option "%s".', option_name);
    end
end
switch lower(unit)
    case 'radian'
        to_doa = @(z) asin(angle(z) / k);
    case 'degree'
        to_doa = @(z) rad2deg(asin(angle(z) / k));
    case 'sin'
        to_doa = @(z) angle(z) / k;
    otherwise
        error('Unkown unit %s', arg_pairs.sampleunit);
end
m = size(R, 1);
if n >= m
    error('Too many sources.');
end
% root MUSIC
[E, l] = eig(R, 'vector');
if isreal(l)
    En = E(:,1:end-n);
else
    % possible complex eigenvalues due to numerical errors.
    [~, idx] = sort(abs(l));
    En = E(:,idx(1:end-n));
end
% compute coeffs
C = En*En';
coeff = zeros(m - 1, 1);
for ii = 1:m-1
    coeff(ii) = sum(diag(C, ii));
end
coeff = [flipud(coeff); sum(diag(C)); conj(coeff)];
% solve
z = roots(coeff); % roots returns a column vector
% find n points inside the unit circle that are also closest to the unit
% circle
% todo: handle the cases when abs(z) == 1
z = z(abs(z) < 1);
[~, idx] = sort(1 - abs(z));

sp = struct();
sp.x_est = sort(to_doa(z(idx(1:n))));
sp.x = sp.x_est;
sp.x_unit = unit;
sp.y = ones(n, 1);
sp.resolved = true;
sp.discrete = true;
end