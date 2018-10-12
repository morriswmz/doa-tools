function design = coprime_1d(mn, d, mode, name)
%COPRIME_1D Generates a 1D co-prime array.
%Syntax:
%   design = COPRIME_1D([2 3], wavelength/2, '2M', 'Co-prime Array');
%Inputs:
%   mn - Co-prime number pair.
%   d - Inter-element spacing.
%   mode - Configuration mode. Can be either 'M', or '2M'. Default is '2M'.
%   name - Custom name of the array. Default is 'Co-prime array (m, n)'.
%Outputs:
%   design - An array design struct.
m = mn(1); n = mn(2);
if m > n
    tmp = m;
    m = n;
    n = tmp;
    warning('M is greater than N. Swapped.');
end
if gcd(m, n) ~= 1
    error('Co-prime number pair expected.');
end
if d <= 0 || ~isreal(d)
    error('d must be a positive real number.');
end
if nargin <= 3
    name = sprintf('Co-prime array (%d, %d)', m, n);
elseif ~ischar(name)
    error('Name must be a string.');
end
if nargin <= 2
    mode = '2m';
end
switch lower(mode)
    case 'm'
        idx = union((0:n-1) * m, (1:m-1) * n);
    case '2m'
        %idx = union((0:n-1) * m, (1:2*m-1) * n);
        idx = [(0:n-1) * m (1:2*m-1) * n];
    otherwise
        error('Unknown mode ''%s''.', mode);
end
design.element_indices = idx;
design.element_positions = design.element_indices*d;
design.element_spacing = d;
design.element_count = length(design.element_indices);
design.dim = 1;
design.type = 'co-prime';
design.name = name;
end