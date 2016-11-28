% array with model errors
clear(); close all;

wavelength = 1; % normalized
d = wavelength / 2;
design_nominal = design_array_1d('coprime', [3 4], d);
% add some position errors
design_perturbed = design_nominal;
design_perturbed.position_errors = 0.05*randn(2, design_perturbed.element_count);
% visualize
visualize_array(design_perturbed, 'VisualizeCoarray', true);