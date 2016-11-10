% using MDL and AIC to find source number, and MUSIC to find DOAs
clear(); close all;

wavelength = 1; % normalized
d = wavelength / 2;
design_ula = design_array_1d('ula', 12, d);
doas = linspace(-pi/3, pi/6, 8);
power_source = 1;
power_noise = 1;
snapshot_count = 100;
source_count = length(doas);

% stochastic (unconditional) model
[~, R] = snapshot_gen_sto(design_ula, doas, wavelength, snapshot_count, power_noise, power_source);

% source number detection
[~, l] = eig(0.5*(R+R'), 'vector');
l = flipud(l);
n_mdl = sn_mdl(l, design_ula.element_count, snapshot_count);
n_aic = sn_aic(l, design_ula.element_count, snapshot_count);
fprintf('There are %d sources.\n', source_count);
fprintf('# of sources estimated by MDL = %d\n', n_mdl);
fprintf('# of sources estimated by AIC = %d\n', n_aic);

% normal MUSIC
sp_normal = music_1d(R, source_count, design_ula, wavelength, 1440);
sp_normal.true_positions = doas;
plot_sp(sp_normal, 'title', 'Normal MUSIC', 'PlotType', 'Polar');

% root MUSIC
sp_root = rmusic_1d(R, source_count, 2*pi*design_ula.element_spacing/wavelength);
sp_root.true_positions = doas;
plot_sp(sp_root, 'Title', 'Root-MUSIC');