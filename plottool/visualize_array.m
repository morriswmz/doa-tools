function hf = visualize_array(design, varargin)
%VISUALIZE_ARRAY Visualizes an array design. If the design struct contains
%'position_errors' field, perturbations will also be visualized.
%Syntax:
%   visualize_array(design);
%Inputs:
%   design - Array design.
%   ... - Options:
%       'ReuseFigure' - Set to true to avoid creating a new figure.
%       'VisualizeCoarray' - Set to true to visualize the difference
%                            coarray.
%       'NominalMarker' - Marker for nominal element positions.
%       'NominalMarkerColor' - Marker color for nominal element positions.
%       'NominalMarkerSize' - Marker size for nominal element positions.
%       'PerturbedMarker' - Marker for for perturbed element positions.
%       'PerturbedMarkerColor' - Marker color for for perturbed element
%                                positions.
%       'PerturbedMarkerSize' - Marker size for perturbed element
%                               positions.
%Output:
%   hf - Handle of the figure.
reuse_figure = false;
visualize_coarray = false;
nominal_marker = 'o';
nominal_marker_color = 'blue';
nominal_marker_size = [];
perturbed_marker = 'x';
perturbed_marker_color = 'red';
perturbed_marker_size = [];
for ii = 1:2:nargin-1
    option_name = varargin{ii};
    option_value = varargin{ii+1};
    switch lower(option_name)
        case 'reusefigure'
            reuse_figure = option_value;
        case 'visualizecoarray'
            visualize_coarray = option_value;
        case 'nominalmarker'
            nominal_marker = option_value;
        case 'nominalmarkercolor'
            nominal_marker_color = option_value;
        case 'nominalmarkersize'
            nominal_marker_size = option_value;
        case 'perturbedmarker'
            perturbed_marker = option_value;
        case 'perturbedmarkercolor'
            perturbed_marker_color = option_value;
        case 'perturbedmarkersize'
            perturbed_marker_size = option_value;
        otherwise
            error('Unknown option name ''%s''.', option_name);
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
visualize_perturbations = isfield(design, 'position_errors');
if visualize_perturbations
    perturbed_positions = get_perturbed_positions(design.element_positions, design.position_errors);
end
% visualize the physical array
plot_array_impl(design.element_positions, ...
    nominal_marker, nominal_marker_color, nominal_marker_size);
if visualize_perturbations
    hold on;
    plot_array_impl(perturbed_positions, ...
        perturbed_marker, perturbed_marker_color, perturbed_marker_size);
    hold off;
    legend('Nominal', 'Perturbed');
end
title(design.name);
axis('equal');
% visualize the coarray
if visualize_coarray
    subplot(2,1,2);
    plot_array_impl(unique_differences(design.element_positions), ...
        nominal_marker, nominal_marker_color, nominal_marker_size);
    if visualize_perturbations
        hold on;
        plot_array_impl(unique_differences(perturbed_positions), ...
            perturbed_marker, perturbed_marker_color, perturbed_marker_size);
        hold off;
        legend('Nominal', 'Perturbed');
    end
    title(['Coarray of ' design.name]);
    axis('equal');
end
end

function plot_array_impl(pos, m, mc, ms)
[dim, n] = size(pos);
if dim == 1
    scatter(pos, zeros(1,n), ms, mc, 'Marker', m, 'MarkerEdgeColor', mc);
elseif dim == 2
    scatter(pos(1,:), pos(2,:), ms, mc, 'Marker', m, 'MarkerEdgeColor', mc);
elseif dim == 3
    scatter3(pos(1,:), pos(2,:), pos(3,:), ms, mc, 'Marker', m, 'MarkerEdgeColor', mc);
else
    error('Incorrect array dimension.');
end
end
