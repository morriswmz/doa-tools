function hf = visualize_array(design, varargin)
%VISUALIZE_ARRAY Visualizes an array design.
%Syntax:
%   visualize_array(design);
%Inputs:
%   design - Array design.
%   ... - Options:
%       'ReuseFigure' - Set to true to avoid creating a new figure.
%       'VisualizeCoarray' - Set to true to visualize the difference
%                            coarray.
%Output:
%   hf - Handle of the figure.
reuse_figure = false;
visualize_coarray = false;
for ii = 1:2:nargin-1
    option_name = varargin{ii};
    option_value = varargin{ii+1};
    switch lower(option_name)
        case 'reusefigure'
            reuse_figure = option_value;
        case 'visualizecoarray'
            visualize_coarray = option_value;
        otherwise
            error('Unknown option name "%s".', option_name);
    end
end
if reuse_figure
    hf = gcf;
else
    hf = figure;
end
if visualize_coarray
    subplot(2,1,1);
end
% visualize the physical array
plot_array_impl(design.element_positions);
title(design.name);
% visualize the coarray
if visualize_coarray
    subplot(2,1,2);
    plot_array_impl(unique_differences(design.element_positions));
    title(['Coarray of ' design.name]);
end
end

function plot_array_impl(pos)
[dim, n] = size(pos);
if dim == 1
    scatter(pos, zeros(1,n));
elseif dim == 2
    scatter(pos(1,:), pos(2,:));
elseif dim == 3
    scatter3(pos(1,:), pos(2,:), pos(3,:));
else
    error('Incorrect array dimension.');
end
axis('equal');
end
