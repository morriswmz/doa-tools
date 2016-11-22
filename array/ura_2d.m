function design = ura_2d(size, d, name)
%URA_2D Generates a 2D uniform rectangular array.
%Syntax:
%   design = URA_2D(size, d[ ,name]);
%Inputs:
%   size - A two element vector determining the size of the array. The
%          resulting array will have size(1) x size(2) elements.
%   d - Inter-element spacing.
%   name - Custom name of the array. Default is 'URA (size(1), size(2))'.
%Outputs:
%   design - An array design struct.
if d <= 0 || ~isreal(d)
    error('d must be a positive real number.');
end
if nargin <= 2
    name = sprintf('URA (%d, %d)', size(1), size(2));
elseif ~ischar(name)
    error('Name must be a string.');
end
design.element_indices = grid2(0, size(1)-1, 0, size(2) - 1, size(1), size(2));
design.element_positions = design.element_indices*d;
design.element_spacing = d;
design.element_count = size(1)*size(2);
design.dim = 2;
design.type = 'ura';
design.name = name;
end

