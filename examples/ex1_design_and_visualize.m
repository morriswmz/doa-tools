% design and visualize arrays
clear(); close all;

wavelength = 1; % normalized
d = wavelength / 2;
design_ula = design_array_1d('ula', 12, d);
design_cp = design_array_1d('coprime', [4 5], d);
design_nested = design_array_1d('nested', [5 7], d);
design_mra = design_array_1d('mra', 12, d);

visualize_array(design_ula, 'VisualizeCoarray', true);
visualize_array(design_cp, 'VisualizeCoarray', true);
visualize_array(design_nested, 'VisualizeCoarray', true);
visualize_array(design_mra, 'VisualizeCoarray', true);