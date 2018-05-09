% Stochastic error model.
% RMSE vs number of snapshots.
% Produces figures similar to Fig. 7 and 8.

%% Configure
clear(); %close all;

wavelength = 1; % normalized
d_0 = wavelength / 2;
design_set = 'SameNumSensor';  % can also be 'SameAperture'
designs = get_design_set(design_set, d_0);
n_designs = length(designs);
if strcmpi(design_set, 'SameNumSensor')
    n_doas = 11;
else
    n_doas = 6;
end
doas = linspace(-pi/3, pi/3, n_doas);

source_power = 1;
noise_power = 1;

perturb_type = 'gaussian';
perturb_std = 0.1*d_0;
n_params_n_snapshots = 20;
params_n_snapshots = floor(linspace(100, 5000, n_params_n_snapshots));
n_repeats = 50; % set to 5000 in the paper

%% Run
mse_em = zeros(n_designs, n_params_n_snapshots, n_repeats);
mse_an = zeros(n_designs, n_params_n_snapshots);

fprintf('Fixed perturbation level, varying number of snapshots:\n');
fprintf('    DOAs: [%s]\n', num2str(doas, '%f '));
fprintf('    SNR: %.1f\n', 10*log10(source_power / noise_power));
fprintf('    # of snapshots: [%d %d]\n', min(params_n_snapshots), max(params_n_snapshots));
fprintf('    Repeat: %d\n', n_repeats);

for dd = 1:n_designs
    cur_design = designs{dd};
    pos_err_info = struct('std', perturb_std, 'mask', ...
        true(cur_design.element_count, 1), 'reference', 'absolute');
    fprintf('Running simulations for %s:\n', cur_design.name);
    progressbar('reset', n_params_n_snapshots);
    for kk = 1:n_params_n_snapshots
        n_snapshots = params_n_snapshots(kk);
        % collect empirical results
        for rr = 1:n_repeats
            [~, R] = snapshot_gen_sto_pos_err(cur_design, doas, wavelength, ...
                n_snapshots, noise_power, source_power, pos_err_info);
            [Rv, ~, ~] = virtual_ula_cov_1d(cur_design, R, 'SS');
            sp = rmusic_1d(Rv, n_doas, 2*pi * cur_design.element_spacing / wavelength);
            mse_em(dd,kk,rr) = sum((sp.x_est - doas).^2) / n_doas;
        end
        % compute analytical results
        c1 = exp(-4 * pi^2 / wavelength^2 * perturb_std^2);
        new_noise_power = (noise_power + (1 - c1) * n_doas * source_power) / c1;
        mse_an(dd,kk) = mean(ecov_coarray_music_1d(cur_design, wavelength, doas, ...
            source_power, new_noise_power, n_snapshots, 'DiagonalsOnly'));
        progressbar('advance');
    end
    progressbar('end');
end

%% Plot
markers = {'x', 'o', 's', 'd', '*', '^'};
hf = figure;
for dd = 1:n_designs
    hp = semilogy(params_n_snapshots, rad2deg(sqrt(mean(squeeze(mse_em(dd,:,:)), 2))), ...
        markers{mod(dd, length(markers))+1}, ...
        'DisplayName', [designs{dd}.name ' empirical']);
    hold on;
    semilogy(params_n_snapshots, rad2deg(sqrt(mse_an(dd,:))), '--',...
        'Color', get(hp, 'Color'), ...
        'DisplayName', [designs{dd}.name ' analytical']);
end
hold off;
grid on;
xlabel('Number of snapshots');
ylabel('RMSE/deg');
legend('show', 'Location', 'northeast');
title('RMSE vs # of Snapshots - Stochastic Error Model');
