function design = ula_1d(n, d, name)
%ULA_1D Generates a 1D ULA.
%Syntax:
%   design = ULA_1D(10, wavelength/2);
%   design = ULA_1D(10, 0.5, 'ULA with 10 sensors');
%Inputs:
%   n - Number of elements.
%   d - Inter-element spacing.
%   name - Custom name of the array. Default is 'ULA with n elements'.
%Outputs:
%   design - An array design struct.
if n <= 0
    error('n_sensor must be a positive integer.');
end
if d <= 0 || ~isreal(d)
    error('d must be a positive real number.');
end
if nargin <= 2
    name = sprintf('ULA with %d elements', n);
elseif ~ischar(name)
    error('Name must be a string.');
end
design.element_indices = (0:n-1);
design.element_positions = design.element_indices*d;
design.element_spacing = d;
design.element_count = n;
design.dim = 1;
design.type = 'ula';
design.name = name;
end

