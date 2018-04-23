% Checks the resolution of SS-MUSIC or DA-MUSIC.
% This script will generate results similar Fig. 2 in the following paper:
% * M. Wang and A. Nehorai, "Coarrays, MUSIC, and the Cramer-Rao Bound,"
%   IEEE Trans. Signal Process., vol. 65, no. 4, pp. 933-946, Feb. 2017.
% This script uses root-MUSIC, will actually yields better resolution
% results than those presented in the original paper. Nevertheless, we
% have similar observations and our analytical results indeeds predicts
% the resolution limit of SS-MUSIC.

clear(); close all;

rv_mode = 'SS';
wavelength = 1;
d_0 = wavelength / 2;
power_source = 1;
power_noise = 1;
n_snapshots = 500;

% designs
design_cp = design_array_1d('coprime', [3 5], d_0, '2M', 'Co-prime Array');
design_nested = design_array_1d('nested', [4 6], d_0, 'Nested Array');
design_mra = design_array_1d('mra', 10, d_0, 'MRA');
designs = { ...
    design_cp, ...
    design_nested, ...
    design_mra
};
n_designs = length(designs);

n_repeats = 500;

theta_0 = pi/6;
delta_theta = linspace(pi/60, pi/600, 20);
n_delta_theta = length(delta_theta);
resolution_success_rate = zeros(n_designs, n_delta_theta);
resolution_threshold_ana = zeros(n_designs);

% run
for dd = 1:n_designs
    design = designs{dd};
    min_delta_ana = delta_theta(end);
    unresolvable = false;
    unresolvable_flags = false(n_delta_theta);
    parfor ii = 1:n_delta_theta
        theta_1 = theta_0 - delta_theta(ii)/2;
        theta_2 = theta_0 + delta_theta(ii)/2;
        doas = [theta_1 theta_2];
        % analytically check the minimal resolvable delta
        unresolvable_flags(ii) = ...
            ~check_resolution_ana(design, wavelength, theta_1, theta_2, ...
                power_source, power_noise, n_snapshots);
        % empirical success rate
        n_success = 0;
        for rr = 1:n_repeats                
            [~, R] = snapshot_gen_sto(design, doas, wavelength, n_snapshots, power_noise, power_source);
            Rv = virtual_ula_cov_1d(design, R, rv_mode);
            sp = rmusic_1d(Rv, 2, 2*pi*d_0/wavelength);
            % requires accurate resolution
            if (check_doa_correctness(doas, sp.x_est, delta_theta(ii)/2))
                n_success = n_success + 1;
            end
        end
        resolution_success_rate(dd, ii) = n_success / n_repeats;
    end
    min_delta_ana_idx = find(unresolvable_flags == true, 1);
    if isempty(min_delta_ana_idx)
        resolution_threshold_ana(dd) = delta_theta(1);
    else
        resolution_threshold_ana(dd) = delta_theta(min_delta_ana_idx);
    end
end

%% plot
figure;
markers = {'+', 'o', 'd', '*'};
res_legends = {};
for dd = 1:n_designs
    hp = plot(rad2deg(delta_theta), resolution_success_rate(dd,:), ['-' markers{dd}]);
    hold on;
    axis([-inf inf 0 1]);
    line(rad2deg([resolution_threshold_ana(dd) resolution_threshold_ana(dd)]), ...
        [0 1], 'Color', hp.Color, 'LineStyle', '--');
    hold on;
    res_legends = [res_legends designs{dd}.name [designs{dd}.name ' Predicted']]; 
end
hold off;
xlabel('\Delta \theta (degree)'); ylabel('Rate of success');
title('Probability of resolution');
legend(res_legends{:});