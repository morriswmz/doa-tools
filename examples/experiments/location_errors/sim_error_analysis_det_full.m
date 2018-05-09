% Deterministic error model.
% Full verification.
% Produces figures similar to Fig. 3.

%% Configure
clear();

wavelength = 1; % normalized
d_0 = wavelength / 2;
designs = get_design_set('SameNumSensor', d_0);

n_designs = length(designs);
n_doas = 11;
doas = linspace(-pi/3, pi/3, n_doas);

source_power = 1;
noise_power = 1;

perturb_type = 'gaussian';
n_params = 20; % 20 is used in the paper.
max_perturb = 0.10*d_0;
max_n_snapshots = 1000;
params_perturb = linspace(0, max_perturb, n_params);
params_n_snapshots = round(linspace(100, max_n_snapshots, n_params));

n_repeats = 1000;
%% Run
mse_em = zeros(n_designs, n_params, n_params, n_repeats);
mse_an = zeros(n_designs, n_params, n_params);

fprintf('Varying position error std and number of snapshots:\n');
fprintf('    DOAs: [%s]\n', num2str(doas, '%f '));
fprintf('    SNR: %.1f\n', 10*log10(source_power / noise_power));
fprintf('    Repeat: %d\n', n_repeats);

for dd = 1:n_designs
    cur_design = designs{dd};
    fprintf('Running simulations for %s:\n', cur_design.name);
    progressbar('reset', n_params^2);
    for ii = 1:n_params
        n_snapshots = params_n_snapshots(ii);
        for kk = 1:n_params
            perturb_std = params_perturb(kk);
            % collect empirical results
            for rr = 1:n_repeats
                pos_err = gen_pos_err(perturb_std, cur_design.element_count, perturb_type);
                perturbed_design = cur_design;
                perturbed_design.position_errors = pos_err;
                [~, R] = snapshot_gen_sto(perturbed_design, doas, wavelength,...
                            n_snapshots, noise_power, source_power);
                [Rv, ~, ~] = virtual_ula_cov_1d(cur_design, R, 'SS');
                sp = rmusic_1d(Rv, n_doas, 2*pi * cur_design.element_spacing / wavelength);
                mse_em(dd,ii,kk,rr) = sum((sp.x_est - doas).^2) / n_doas;
            end
            % compute analytical results
            [~, perturb_cov] = gen_pos_err(perturb_std, cur_design.element_count, perturb_type);
            error_stat = struct;
            error_stat.PositionErrorCov = perturb_cov;
            mse_an(dd,ii,kk) = sum(ecov_perturbed_coarray_music_1d(cur_design, wavelength, ...
                doas, source_power, noise_power, n_snapshots, error_stat, 'DiagonalsOnly')) / n_doas;
            progressbar('advance');
        end
    end
    progressbar('end');
end

%% Plot
rmse_em = rad2deg(sqrt(mean(mse_em, 4)));
rmse_an = rad2deg(sqrt(mse_an));
rel_err = abs(rmse_em - rmse_an) ./ rmse_an;

for dd = 1:n_designs
    figure;
    imagesc(params_perturb/d_0, params_n_snapshots, squeeze(rel_err(dd,:,:))); axis ij;
    xlabel('\delta/d_0');
    ylabel('Number of Snapshots');
    title(designs{dd}.name);
    colormap gray;
    colormap(flipud(colormap));
    colorbar;
    set(gca, 'Clim', [0 0.5]);
    set(gca,'YDir','normal');
end
