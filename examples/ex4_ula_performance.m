% compares the performance of ULA with the CRB.
clear(); close all;

wavelength = 1; % normalized
d = wavelength / 2;
design_ula = design_array_1d('ula', 12, d);
doas = linspace(-pi/3, pi/4, 8);
power_source = 1;
snapshot_count = 200;
source_count = length(doas);
P = eye(source_count)*power_source; % equally powered sources

param_power_noise = 10.^linspace(-2, 2, 20);
n_param = length(param_power_noise);
n_repeat = 200; % 200 Monte Carlo runs for each parameter
mse = zeros(n_param, 1);
crb_sto = zeros(n_param, 1);
crb_det = zeros(n_param, 1);
crb_sto_uc = zeros(n_param, 1);
snr = zeros(n_param, 1);
progressbar('reset', n_param);
for ii = 1:n_param
    power_noise = param_power_noise(ii);
    snr(ii) = 10*log10(power_source / power_noise);
    cur_se = 0;
    cur_crb_det = 0;
    for rr = 1:n_repeat
        [~, R, S] = snapshot_gen_sto(design_ula, doas, wavelength, snapshot_count, power_noise, power_source);
        sp = rmusic_1d(R, source_count, 2*pi*design_ula.element_spacing/wavelength);
        cur_se = cur_se + sum((sp.x_est - doas).^2);
        P_est = (S*S')/snapshot_count;
        cur_crb_det = cur_crb_det + mean(diag(...
            crb_general_det_1d(design_ula, wavelength, doas, P_est, power_noise, snapshot_count)));
    end
    mse(ii) = cur_se / (source_count * n_repeat);
    crb_sto(ii) = mean(diag(crb_general_sto_1d(design_ula, wavelength, doas, P, power_noise, snapshot_count)));
    crb_det(ii) = cur_crb_det / n_repeat;
    crb_sto_uc(ii) = mean(diag(crb_uc_sto_1d(design_ula, wavelength, doas, power_source, power_noise, snapshot_count)));
    progressbar('advance');
end
fprintf('\n');
% plot
% we should observe that the MSE approaches the CRB in high SNR regions,
% and that the stochastic CRB is tighter than the deterministic CRB.
% With the additional assumption of uncorrelated sources, we expect a lower
% CRB.
% As SNR -> infinity, all three CRBs converges to the same one.
semilogy(snr, mse, '-x', snr, crb_sto, '--', snr, crb_det, '--', ...
    snr, crb_sto_uc, '--');
xlabel('SNR (dB)'); ylabel('MSE / rad^2'); grid on;
legend('MSE', 'Stochastic CRB', 'Deterministic CRB', 'Stochastic CRB (Uncorrelated)');
title('MSE vs. CRB');
