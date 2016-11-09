% resolve more uncorrelated sources than the number of array elements using
% co-prime arrays
clear(); close all;

wavelength = 1; % normalized
d = wavelength / 2;
% this co-prime array only have 9 elements
design_cp = design_array_1d('coprime', [3 4], d);
% but we have 10 sources here
doas = linspace(-pi/3, pi/3, 10);
power_source = 1;
power_noise = 1;
snapshot_count = 100;
source_count = length(doas);

% unconditional model
A = steering_matrix_1d(design_cp, wavelength, doas);
X = sqrt(power_source) * randccsn(source_count, snapshot_count);
N = sqrt(power_noise) * randccsn(design_cp.element_count, snapshot_count);
Y = A*X + N;
R = (Y*Y') / snapshot_count;

% SS-MUSIC
[Rss, dss] = virtual_ula_cov_1d(design_cp, R, 'SS');
sp = music_1d(Rss, source_count, dss, wavelength, 1440);
sp.true_positions = doas;
plot_sp(sp, 'title', ['SS-MUSIC using ' design_cp.name]);
