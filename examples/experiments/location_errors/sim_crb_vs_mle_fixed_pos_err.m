% Demonstrates the CRB obtained from (33) is indeed attainable by the stochastic
% MLE.
% Produces figures similar to Fig. 9.

%% Configure
clear();
rng(42);

wavelength = 1; % normalized
d_0 = wavelength / 2;
design_set = 'SameNumSensor'; % must have the same number of sensors!
designs = get_design_set(design_set, d_0);
n_designs = length(designs);
n_snapshots = 5000;
n_doas = 11;
% Note: here these numbers are reduces because the MLE takes a long time to
% compute. Consider increasing them to obtain a more stable result.
n_repeats = 10;
n_params = 10;
doas = linspace(-pi/3, pi/3, n_doas);

power_source = 1;
params_SNR = linspace(-20, 20, n_params);

% generate fixed position error vectors
pos_err = randn(2, designs{1}.element_count) * (0.1 * d_0);
pos_err(:,[1 2]) = 0;
mask_u = pos_err(1,:) ~= 0;
mask_v = pos_err(2,:) ~= 0;
pos_err_mask = [mask_u; mask_v];

%% Run
mses = zeros(n_designs, n_params, n_repeats);
crbs = zeros(n_designs, n_params);
crbs_nominal = zeros(n_designs, n_params);

for dd = 1:n_designs
    design_nominal = designs{dd};
    design = design_nominal;
    design.position_errors = pos_err(:,1:design.element_count);
    fprintf('Running simulation for ''%s''\n', design.name);
    
    % Remember to remove the progressbar calls when using parfor.
    progressbar('reset', n_params * n_repeats);
    for ii = 1:n_params
        SNR = params_SNR(ii);
        power_noise = power_source * 10^(-SNR/10);
        for rr = 1:n_repeats
            [~, R] = snapshot_gen_sto(design, doas, wavelength, n_snapshots, power_noise, power_source);
            [doa_est, p_est, u_est, v_est, noise_est] = mle_uc_sto_pos_err(R, n_doas, design, wavelength, [], pos_err_mask);
            mses(dd, ii, rr) = mean((doa_est(:) - doas(:)).^2);
            progressbar('advance');
        end
        CRB = crb_uc_sto_pos_err(design, wavelength, doas, power_source, power_noise, n_snapshots, pos_err_mask);
        crbs(dd, ii) = mean(diag(CRB));
        CRB_nominal = crb_uc_sto_1d(design_nominal, wavelength, doas, power_source, power_noise, n_snapshots);
        crbs_nominal(dd, ii) = mean(diag(CRB_nominal));
    end
    progressbar('end');
end
    
%% Plot
% The MLE obtained from fmincon may produce very large errors, leading to
% outliers which will greatly affect the empirical MSE. Set the following
% variable to true to suppress DOA estimates that are too large.
suppress_outliers = false;
processed_mses = mses;
if suppress_outliers
    max_mse = ((max(doas) - min(doas)) / (n_doas - 1) / 2)^2;
    processed_mses(processed_mses > max_mse) = nan;
end
mses_final = nanmean(processed_mses, 3);

for dd = 1:n_designs
    hf = figure;
    semilogy(params_SNR, mses_final(dd,:), 'kx', ...
        params_SNR, crbs(dd,:), 'k--', params_SNR, crbs_nominal(dd,:), 'k-.');
    legend('MSE', 'CRB', 'CRB (location-error free)');
    axis([-inf inf -inf 1]);
    xlabel('SNR (dB)'); ylabel('MSE (rad^2)');
    title(designs{dd}.name);
    fig_pos = get(hf, 'Position');
    set(hf, 'Position', [fig_pos(1) fig_pos(2) 400 300]);
    set(hf, 'Color', 'White');
end
