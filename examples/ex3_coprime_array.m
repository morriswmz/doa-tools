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

% stochastic (unconditional) model
[~, R] = snapshot_gen_sto(design_cp, doas, wavelength, snapshot_count, power_noise, power_source);

% SS-MUSIC and DA-MUSIC
[Rss, dss] = virtual_ula_cov_1d(design_cp, R, 'SS');
[Rda, dda] = virtual_ula_cov_1d(design_cp, R, 'DA');
sp_ss = music_1d(Rss, source_count, dss, wavelength, 1440);
sp_da = music_1d(Rda, source_count, dda, wavelength, 1440);
sp_ss.true_positions = doas;
sp_da.true_positions = doas;
figure;
subplot(2,1,1);
plot_sp(sp_ss, 'Title', ['SS-MUSIC using ' design_cp.name], 'ReuseFigure', true);
subplot(2,1,2);
plot_sp(sp_da, 'Title', ['DA-MUSIC using ' design_cp.name], 'ReuseFigure', true);

% apply MVDR directly
sp_mvdr = mvdr_1d(R, source_count, design_cp, wavelength, 1440);
sp_mvdr.true_positions = doas;
plot_sp(sp_mvdr, 'Title', ['Direct MVDR using ' design_cp.name]);