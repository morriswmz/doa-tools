% Deterministic error model.
% MSE vs. SNR under different levels of perturbations.
% Produces figures similar to Fig. 6.

%% Configure
clear();

wavelength = 1; % normalized
d_0 = wavelength / 2;
design = design_array_1d('coprime', [2 3], d_0, '2M', 'Co-prime (2,3)');

n_doas = 6;
doas = linspace(-pi/3, pi/3, n_doas);

n_snapshots = 5000;
source_power = 1;
noise_power = 1;

perturb_type = 'gaussian';
params_perturb = [0 0.01 0.02 0.04] * d_0;
n_params_perturb = length(params_perturb);
max_perturb = max(params_perturb);
n_params_snr = 20;
params_snr = linspace(-15, 15, n_params_snr);
n_repeat = 1000;

%% Run

mse_em_ss = zeros(n_params_perturb, n_params_snr, n_repeat);
mse_an = zeros(n_params_perturb, n_params_snr);

fprintf('Fixed snapshot count, varying SNR parameters:\n');
fprintf('    DOAs: [%s]\n', num2str(doas, '%f '));
fprintf('    Number of snapshots: %d\n', n_snapshots);
fprintf('    Perturbation range: [0 %f]d_0\n', max_perturb/d_0);
fprintf('    Repeat: %d\n', n_repeat);

progressbar('reset', n_params_perturb * n_params_snr);
for kk = 1:n_params_perturb
    perturb_std = params_perturb(kk);
    for ss = 1:n_params_snr
        noise_power = 10^(-params_snr(ss) / 10);
        % collect empirical results
        for rr = 1:n_repeat
            pos_err = gen_pos_err(perturb_std, design.element_count, perturb_type);
            perturbed_design = design;
            perturbed_design.position_errors = pos_err;
            [~, R] = snapshot_gen_sto(perturbed_design, doas, wavelength, n_snapshots, noise_power, source_power);
            [Rv, ~, ~] = virtual_ula_cov_1d(design, R, 'SS');
            sp = rmusic_1d(Rv, n_doas, 2*pi * design.element_spacing / wavelength);
            mse_em_ss(kk, ss, rr) = sum((sp.x_est - doas).^2) / n_doas;
        end
        % compute analytical results
        [~, perturb_cov] = gen_pos_err(perturb_std, design.element_count, perturb_type);
        error_stat = struct;
        error_stat.PositionErrorCov = perturb_cov;
        mse_an(kk, ss) = sum(ecov_perturbed_coarray_music_1d(design, wavelength, ...
            doas, source_power, noise_power, n_snapshots, error_stat, 'DiagonalsOnly')) / n_doas;
        progressbar('advance');
    end
end
progressbar('end');

%% Plot
markers = {'x', 'o', 's', 'd', '*', '^'};
figure;
for pp = 1:n_params_perturb
    perturb_str = sprintf('\\sigma_p/d_0 = %.2f', params_perturb(pp)/d_0);
    hp = semilogy(params_snr, rad2deg(sqrt(nanmean(squeeze(mse_em_ss(pp,:,:)), 2))), ...
        ['' markers{mod(pp, length(markers))+1}], ...
        'DisplayName', ['SSM empirical:' perturb_str]);
    hold on;
    semilogy(params_snr, rad2deg(sqrt(mse_an(pp,:))), '--',...
        'Color', get(hp, 'Color'), ...
        'DisplayName', ['SSM analytical:' perturb_str]);
    hold on;
end
xlabel('SNR/dB');
ylabel('RMSE/deg');
hold off;
legend('show');
title('RMSE vs SNR - Deterministic Error Model');
