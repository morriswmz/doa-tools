function sp = rmusic_1d(R, n, k, varargin)
%RMUSIC_1D 1D root MUSIC for ULAs.
%Syntax:
%   sp = RMUSIC_1D(R, n, k);
%Inputs:
%   R - Sample covariance matrix of a ULA.
%   n - Number of sources.
%   k - 2*pi*inter_element_spacing/wavelength
%   ... - Options:
%       'Unit' - Can be 'radian', 'degree', or 'sin'. Default value is
%                'radian'.
%Output:
%   sp - Spectrum structure with the following fields:
%           x - An 1 x grid_size vector.
%           y - An 1 x grid_size vector. Calling `plot(x, y)` will plot the
%               spectrum.
%           x_est - An 1 x n vector storing the estimated DOAs.
%           x_unit - The same as the unit specified by 'Unit'.
%           resolved - Constant value true.
%           discrete - Constant value true.
%Reference:
%   [1] H. L. Van Trees, Optimum array processing. New York: Wiley, 2002.
unit = 'radian';
for ii = 1:2:nargin-3
    option_name = varargin{ii};
    option_value = varargin{ii+1};
    switch lower(option_name)
        case 'unit'
            unit = option_value;
        otherwise
            error('Unknow option ''%s''.', option_name);
    end
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
z = roots(coeff)'; % roots returns a column vector
% find n points inside the unit circle that are also closest to the unit
% circle
nz = length(z);
mask = true(nz, 1);
for ii = 1:nz
    absz = abs(z(ii));
    if absz > 1
        mask(ii) = false;
    elseif absz == 1
        % find the closest point and remove it
        idx = -1;
        dist = inf;
        for jj = 1:nz
            if jj ~= ii && mask(jj)
                cur_dist = abs(z(ii) - z(jj));
                if cur_dist < dist
                    dist = cur_dist;
                    idx = jj;
                end
            end
        end
        mask(idx) = false;
    end
end
z = z(mask);
[~, idx] = sort(1 - abs(z));

sp = struct();
sp.x_est = sort(-cm2doa(z(idx(1:n)), k, unit));
sp.x = sp.x_est;
sp.x_unit = unit;
sp.y = ones(1, n);
sp.resolved = true;
sp.discrete = true;
end