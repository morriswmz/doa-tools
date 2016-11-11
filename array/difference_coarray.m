function dco = difference_coarray(design, name)
%DIFFERENCE_COARRAY Obtains the full difference coarray of the given array.
%Syntax:
%   dco = difference_coarray(design, 'name');
%Input:
%   design - Array design.
%   name - (Optional) name of the array.
%Output:
%   dco - Difference coarray.
if nargin <= 1
    name = 'Difference coarray';
end
if isfield(design, 'element_spacing') && isfield(design, 'element_indices')
    dco.element_indices = unique_differences(design.element_indices);
    dco.element_spacing = design.element_spacing;
    dco.element_positions = dco.element_indices * design.element_spacing;
else
    dco.element_positions = unique_differences(design.element_positions);
end
dco.element_count = size(dco.element_positions, 2);
dco.type = 'custom';
dco.dim = design.dim;
dco.name = name; 
end

