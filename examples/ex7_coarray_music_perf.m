% performance evaluation of coarray-based MUSIC
clear(); close all;

% set up
wavelength = 1; % normalized
d = wavelength / 2;
design = design_array_1d('coprime', [3 4], d);
doas = linspace(-pi/3, pi/3, 10);
power_source = 1;
snapshot_count = 200;
source_count = length(doas);

% we vary SNR from -20dB to 20dB
n_param_power_noise = 20;
param_power_noise = 10.^(-linspace(-20, 20, n_param_power_noise)/10);
n_repeat = 200;
mse_em = zeros(n_param_power_noise, 1);
mse_an = zeros(n_param_power_noise, 1);
crb = zeros(n_param_power_noise, 1);
snr = zeros(n_param_power_noise, 1);
progressbar('reset', n_param_power_noise);
for ii = 1:n_param_power_noise
    power_noise = param_power_noise(ii);
    snr(ii) = 10*log10(power_source / power_noise);
    cur_se = 0;
    for rr = 1:n_repeat
        [~, R, ~] = snapshot_gen_sto(design, doas, wavelength, snapshot_count, power_noise, power_source);
        [Rss, dss] = virtual_ula_cov_1d(design, R, 'SS');
        sp = rmusic_1d(Rss, source_count, 2*pi*design.element_spacing/wavelength);
        cur_se = cur_se + sum((sp.x_est - doas).^2);
    end
    mse_em(ii) = cur_se / (source_count * n_repeat);
    mse_an(ii) = mean(ecov_coarray_music_1d(design, wavelength, doas, power_source, power_noise, snapshot_count, 'DiagonalsOnly'));
    crb(ii) = mean(diag(crb_uc_sto_1d(design, wavelength, doas, power_source, power_noise, snapshot_count)));
    progressbar('advance');
end
fprintf('\n');

% plot
semilogy(snr, mse_em, '-x', snr, mse_an, '-', snr, crb, '--');
xlabel('SNR (dB)'); ylabel('MSE / rad^2'); grid on;
legend('Empirical MSE', 'Analytical MSE', 'Stochastic CRB (Uncorrelated)');
title('MSE of SS-MUSIC vs. CRB');