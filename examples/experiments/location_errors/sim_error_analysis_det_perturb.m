% Deterministic error model.
% Demonstrates how the MSE various as the perturbation level increases.
% Produces figures similar to Fig. 4, 5.

%% Configure
clear();

wavelength = 1; % normalized
d_0 = wavelength / 2;
design_set = 'SameAperture'; % can also be 'SameAperture'
designs = get_design_set(design_set, d_0);

n_designs = length(designs);
if strcmpi(design_set, 'SameNumSensor')
    n_doas = 11;
else
    n_doas = 6;
end
doas = linspace(-pi/3, pi/3, n_doas);

n_snapshots = 1000;
source_power = 1;
noise_power = 1;

perturb_type = 'gaussian';
n_params_perturb = 20;
max_perturb = 0.08*d_0;
params_perturb = linspace(0, max_perturb, n_params_perturb);
n_repeats = 1000;

%% Run
mse_em = zeros(n_designs, n_params_perturb, n_repeats);
mse_an = zeros(n_designs, n_params_perturb);

fprintf('Fixed snapshot count, varying pos err std parameters:\n');
fprintf('    DOAs: [%s]\n', num2str(doas, '%f '));
fprintf('    Number of snapshots: %d\n', n_snapshots);
fprintf('    SNR: %.1f\n', 10*log10(source_power / noise_power));
fprintf('    Perturbation range: [0 %f]d_0\n', max_perturb/d_0);
fprintf('    Repeat: %d\n', n_repeats);

for dd = 1:n_designs
    cur_design = designs{dd};
    fprintf('Running simulations for %s:\n', cur_design.name);
    progressbar('reset', n_params_perturb);
    for kk = 1:n_params_perturb
        % collect empirical results
        cur_mses = zeros(n_repeats, 1);
        perturb_std = params_perturb(kk);
        for rr = 1:n_repeats
            pos_err = gen_pos_err(perturb_std, cur_design.element_count, perturb_type);
            perturbed_design = cur_design;
            perturbed_design.position_errors = pos_err;
            [~, R] = snapshot_gen_sto(perturbed_design, doas, wavelength, n_snapshots, noise_power, source_power);
            [Rv, ~, ~] = virtual_ula_cov_1d(cur_design, R, 'SS');
            sp = rmusic_1d(Rv, n_doas, 2*pi * cur_design.element_spacing / wavelength);
            mse_em(dd,kk,rr) = sum((sp.x_est - doas).^2) / n_doas;
        end
        % compute analytical results
        [~, perturb_cov] = gen_pos_err(params_perturb(kk), cur_design.element_count, perturb_type);
        error_stat = struct;
        error_stat.PositionErrorCov = perturb_cov;
        mse_an(dd, kk) = sum(ecov_perturbed_coarray_music_1d(cur_design, wavelength, ...
            doas, source_power, noise_power, n_snapshots, error_stat, 'DiagonalsOnly')) / n_doas;
        progressbar('advance');
    end
    progressbar('end');
end

%% Plot
markers = {'x', 'o', 's', 'd', '*', '^', '>', '<', 'h'};
figure;
for dd = 1:n_designs
    hp = plot(params_perturb/d_0, rad2deg(sqrt(mean(squeeze(mse_em(dd,:,:)), 2))), ...
        ['' markers{mod(dd, length(markers))+1}], ...
        'DisplayName', [designs{dd}.name ' empirical']);
    hold on;
    plot(params_perturb/d_0, rad2deg(sqrt(mse_an(dd,:))), '--', 'Color', get(hp, 'Color'), ...
        'DisplayName', [designs{dd}.name ' analytical']);
end
hold off;
xlabel('\delta_p/d_0');
ylabel('RMSE/deg');
grid on;
legend('show');
legend('Location', 'northwest');
title('RMSE vs Perturbation Level - Deterministic Error Model');
