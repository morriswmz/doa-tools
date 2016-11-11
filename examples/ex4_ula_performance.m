% compares the performance of ULA with the CRB.
clear(); close all;

wavelength = 1; % normalized
d = wavelength / 2;
design_ula = design_array_1d('ula', 12, d);
doas = linspace(-pi/3, pi/4, 8);
power_source = 1;
snapshot_count = 100;
source_count = length(doas);
P = eye(source_count)*power_source;

param_power_noise = 10.^linspace(-2, 2, 20);
n_param = length(param_power_noise);
n_repeat = 200; % 200 Monte Carlo runs for each parameter
mse = zeros(n_param, 1);
crb = zeros(n_param, 1);
snr = zeros(n_param, 1);
progressbar('reset', n_param);
for ii = 1:n_param
    power_noise = param_power_noise(ii);
    snr(ii) = 10*log10(power_source / power_noise);
    cur_se = 0;
    for rr = 1:n_repeat
        [~, R] = snapshot_gen_sto(design_ula, doas, wavelength, snapshot_count, power_noise, power_source);
        sp = rmusic_1d(R, source_count, 2*pi*design_ula.element_spacing/wavelength);
        cur_se = cur_se + sum((sp.x_est - doas).^2);
    end
    mse(ii) = cur_se / (source_count * n_repeat);
    crb(ii) = mean(diag(crb_general_sto_1d(design_ula, wavelength, doas, P, power_noise, snapshot_count)));
    progressbar('advance');
end
fprintf('\n');
% plot
semilogy(snr, mse, '-x', snr, crb, '--');
xlabel('SNR (dB)'); ylabel('MSE / rad^2');
legend('MSE', 'CRB');
title('MSE vs. CRB');
