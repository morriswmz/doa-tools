% array with model errors
clear(); close all;

wavelength = 1; % normalized
d = wavelength / 2;
design_nominal = design_array_1d('ula', 10, d);
% add some gain errors
design_perturbed = design_nominal;
design_perturbed.gain_errors = 1.0 + sqrt(0.1)*randn(design_perturbed.element_count, 1);
% visualize
visualize_array(design_perturbed, 'VisualizeCoarray', true);

% stochastic model
doas = linspace(-pi/3, pi/3, 5);
power_source = 1;
power_noise = 1;
snapshot_count = 1000;
source_count = length(doas);

[~, R] = snapshot_gen_sto(design_perturbed, doas, wavelength, snapshot_count, ...
                            power_noise, power_source);
% use MUSIC with nominal steering vector
sp_unknown_perturb = music_1d(R, source_count, design_nominal, wavelength, 180, ...
                    'RefineEstimates', true);
sp_unknown_perturb.real_positions = doas;
% use MUSIC with known perturbations steering vector
sp_known_perturb = music_1d(R, source_count, design_perturbed, wavelength, 180, ...
                    'RefineEstimates', true);
sp_known_perturb.real_positions = doas;

figure;
subplot(2,1,1);
plot_sp(sp_unknown_perturb, 'ReuseFigure', true);
subplot(2,1,2);
plot_sp(sp_known_perturb, 'ReuseFigure', true);

rmse = @(real, est) sqrt(mean((real - est).^2));
fprintf('RMSE perturbations unknown = %.3e\n', rmse(doas, sp_unknown_perturb.x_est));
fprintf('RMSE perturbations known   = %.3e\n', rmse(doas, sp_known_perturb.x_est));