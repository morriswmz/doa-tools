function design = nested_1d(mn, d, name)
%NESTED_1D Generates a 1D nested array.
%Syntax:
%   design = NESTED_1D([2 3], wavelength/2, 'Nested Array');
%Inputs:
%   mn - Parameter pair.
%   d - Inter-element spacing.
%   name - Custom name of the array. Default is 'Nested array (m, n)'.
%Outputs:
%   design - An array design struct.
n1 = mn(1); n2 = mn(2);
if d <= 0 || ~isreal(d)
    error('d must be a positive real number.');
end
if nargin <= 2
    name = sprintf('Nested array (%d, %d)', n1, n2);
elseif ~ischar(name)
    error('Name must be a string.');
end
design.element_indices = union(1:n1, (n1+1)*(1:n2)) - 1;
design.element_positions = design.element_indices*d;
design.element_spacing = d;
design.element_count = length(design.element_indices);
design.dim = 1;
design.type = 'nested';
design.name = name;
end