% checks the resolution limit of root-MUSIC
clear(); close all;

wavelength = 1; % normalized
d = wavelength / 2;
design = design_array_1d('ula', 10, d); % 10-element ULA
power_source = 1;
power_noise = 1;
snapshot_count = 100;

delta_theta = deg2rad(linspace(0.05, 5, 20)); % 0.1 deg - 10 deg
n_delta_theta = length(delta_theta);
n_repeat = 500;
res_prob = zeros(n_delta_theta, 1);
progressbar('reset', n_delta_theta);
for ii = 1:n_delta_theta
    n_resolved = 0;
    for rr = 1:n_repeat
        doas = [-delta_theta(ii)/2 delta_theta(ii)/2];
        [~, R] = snapshot_gen_sto(design, doas, wavelength, snapshot_count, power_noise, power_source);
        sp = rmusic_1d(R, 2, 2*pi*design.element_spacing/wavelength);
        % resolution condition:
        %   theta_1 \in (-delta_theta, 0)
        %   theta_2 \in (0, delta_theta)
        % This condition is pretty strict. It requires the estimated DOAs
        % are reasonably close to the correct ones.
        % You may change it to a milder condition.
        if sp.x_est(1) < 0 && sp.x_est(1) > -delta_theta(ii) ...
                && sp.x_est(2) > 0 && sp.x_est(2) < delta_theta(ii)
            n_resolved = n_resolved + 1;
        end
    end
    res_prob(ii) = n_resolved / n_repeat;
    progressbar('advance');
end
fprintf('\n');

plot(rad2deg(delta_theta), res_prob); grid on;
xlabel('Source separation (deg)'); ylabel('Resolution Probability');