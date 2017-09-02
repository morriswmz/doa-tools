function sp = esprit_1d(R, n, k, varargin)
%ESPRIT_1D 1D ESPRIT for ULAs.
%Syntax:
%   sp = ESPRIT_1D(R, n, design, wavelength, grid_size, ...);
%Inputs:
%   R - Sample covariance matrix.
%   n - Number of sources.
%   k - 2*pi*inter_element_spacing/wavelength.
%   ... - Options:
%           'Unit' - Can be 'radian', 'degree', or 'sin'. Default value is
%                   'radian'.
%           'Displacement' - The displacement between the two overlapping
%                            subarrays measured in number of element
%                            spacings. Default value is 1.
%                            Note: increasing this value will lead to
%                            smaller unambiguous range. Make sure your DOAs
%                            falls within the unambiguous range.
%           'Formulation' - Either 'TLS' (Total Lease Squares) or 'LS'
%                           (Least Squares). Default value is 'TLS'.
%           'RowWeights' - Specifies the row weights with a vector or
%                          string. Default value is 'Default', which
%                          generates the following weight vector:
%                           [1 sqrt(2) sqrt(3) ... sqrt(3) sqrt(2) 1]
%                          You can disable row weighting by passing in
%                          'Identity' or 'Off'.
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
use_tls = true;
ds = 1;
row_weights = [];
use_row_weights = false;
unit = 'radian';
for ii = 1:2:nargin-3
    option_name = varargin{ii};
    option_value = varargin{ii+1};
    switch lower(option_name)
        case 'unit'
            unit = option_value;
        case 'formulation'
            switch lower(option_value)
                case 'ls'
                    use_tls = false;
                case 'tls'
                    use_tls = true;
                otherwise
                    error('Formulation must be either ''LS'' or ''TLS''.');
            end
        case 'displacement'
            if option_value < 1 || mod(option_value, 1) ~= 0
                error('Displacement must be an integer that is greater or equal to one.');
            end
            ds = option_value;
        case 'rowweights'
            if ischar(option_value)
                switch lower(option_value)
                    case 'default'
                        use_row_weights = true;
                    case {'off', 'identity'}
                        use_row_weights = false;
                    otherwise
                        error('Either specify the row weights manually or pass in ''Default'' to use the default weights, or ''Identity'', ''Off'' to disable row weighting.');
                end
            else
                use_row_weights = true;
                row_weights = option_value(:);
            end
        otherwise
            error('Unknow option "%s".', option_name);
    end
end
m = size(R, 1);
if n > m - ds
    error('Too many sources.');
end
if ~isempty(row_weights) && length(row_weights) ~= m - ds
    error('The dimension of the row weights vector is not equal to (m - displacement).');
end
% ESPRIT
[E, ~] = eig(0.5*(R + R'), 'vector');
Es = E(:,end - n + 1:end);
% apply weights if necessary
if use_row_weights
    if isempty(row_weights)
        % default weights
        if mod(m - ds, 2) == 1
            w_max = (m - ds - 1)/2;
            row_weights = sqrt([1:w_max w_max + 1 w_max:-1:1])';
        else
            w_max = (m - ds)/2;
            row_weights = sqrt([1:w_max w_max:-1:1])';
        end
    end
    Es1 = bsxfun(@times, row_weights, Es(1:end - ds,:));
    Es2 = bsxfun(@times, row_weights, Es(ds + 1:end,:));
else
    Es1 = Es(1:end - ds,:);
    Es2 = Es(ds + 1:end,:);
end
if use_tls
    % TLS estimate
    C = [Es1 Es2];
    C = C'*C;
    C = 0.5*(C + C');
    [V, l] = eig(C, 'vector');
    [~, idx] = sort(real(l), 'descend');
    V = V(:,idx);
    V12 = V(1:n,n + 1:end);
    V22 = V(n + 1:end,n + 1:end);
    Phi = -V12/V22;
else
    % LS estimate
    Phi = (Es1'*Es1)\(Es1'*Es2);
end
% convert z to spectrum
z = eig(Phi);
sp = struct();
sp.x_est = sort(cm2doa(z, k*ds, unit));
sp.x = sp.x_est;
sp.x_unit = unit;
sp.y = ones(1, n);
sp.resolved = true;
sp.discrete = true;
end

