% design and visualize arrays
clear(); close all;

wavelength = 1; % normalized
d = wavelength / 2;
design_ula = design_array_1d('ula', 8, d);
design_cp = design_array_1d('coprime', [3, 4], d);

visualize_array(design_ula, 'VisualizeCoarray', true);
visualize_array(design_cp, 'VisualizeCoarray', true);