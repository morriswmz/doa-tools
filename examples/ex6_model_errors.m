% array with model errors
clear(); close all;

wavelength = 1; % normalized
d = wavelength / 2;
n_sensor = 10;
design_nominal = design_array_1d('ula', n_sensor, d);
% add some gain errors
design_perturbed = design_nominal;
design_perturbed.gain_errors = 1.0 + sqrt(0.1)*randn(n_sensor, 1);
% add some phase errors
design_perturbed.phase_errors = exp(1j*sqrt(0.1)*randn(n_sensor, 1));
% add some position errors
pos_err = randn(2, n_sensor, 1) * sqrt(0.01*d);
pos_err(1,2:end) = pos_err(1,2:end) - pos_err(1,1);
pos_err(2,2:end) = pos_err(2,2:end) - pos_err(2,1);
pos_err(:,1) = 0;
design_perturbed.position_errors = pos_err;
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
% MUSIC with nominal steering vector (assuming no perturbations)
sp_unknown_perturb = music_1d(R, source_count, design_nominal, wavelength, 180, ...
                    'RefineEstimates', true);
sp_unknown_perturb.true_positions = doas;
% MUSIC with known perturbations
sp_known_perturb = music_1d(R, source_count, design_perturbed, wavelength, 180, ...
                    'RefineEstimates', true);
sp_known_perturb.true_positions = doas;

figure;
subplot(2,1,1);
plot_sp(sp_unknown_perturb, 'ReuseFigure', true, 'Title', 'Perturbations Unknown');
subplot(2,1,2);
plot_sp(sp_known_perturb, 'ReuseFigure', true, 'Title', 'Perturbations Known');

rmse = @(real, est) sqrt(mean((real - est).^2));
fprintf('RMSE with unknown perturbations = %.3e\n', rmse(doas, sp_unknown_perturb.x_est));
fprintf('RMSE with known perturbations   = %.3e\n', rmse(doas, sp_known_perturb.x_est));