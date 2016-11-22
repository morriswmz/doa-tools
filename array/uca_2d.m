function design = uca_2d(n, r, name)
%UCA_2D Generates a uniform circular array centered at the origin.
%Syntax:
%   design = UCA_2D(n, r[, name]);
%Inputs:
%   n - Number of elements.
%   r - Radius.
%   name - Custom name of the array. Default is 'UCA with n elements'.
%Outputs:
%   design - An array design struct.
if n <= 0
    error('n_sensor must be a positive integer.');
end
if nargin <= 2
    name = sprintf('UCA with %d elements', n);
elseif ~ischar(name)
    error('Name must be a string.');
end
step = 2*pi/n;
theta = 0:step:(2*pi-step);
design.element_positions = r*[cos(theta); sin(theta)];
design.element_spacing = 2*r*sin(step);
design.element_count = n;
design.dim = 2;
design.type = 'uca';
design.name = name;
end

