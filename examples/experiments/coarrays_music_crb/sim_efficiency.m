% Efficiency analysis for SS-MUSIC (or DA-MUSIC).
% This script will generate Fig. 4 in the following paper:
% * M. Wang and A. Nehorai, "Coarrays, MUSIC, and the Cramer-Rao Bound,"
%   IEEE Trans. Signal Process., vol. 65, no. 4, pp. 933-946, Feb. 2017.

clear(); close all;

n_snapshots = 1; % normalized
wavelength = 1;
d_0 = wavelength / 2;
power_source = 1;

design_cp = design_array_1d('coprime', [3 5], d_0, '2M', 'Co-prime Array');
design_nested = design_array_1d('nested', [4 6], d_0, 'Nested Array');
design_mra = design_array_1d('mra', 10, d_0, 'MRA');
designs = { ...
    design_cp, ...
    design_nested, ...
    design_mra
};
n_designs = length(designs);

SNRs = linspace(-10, 20, 30);
n_SNRs = length(SNRs);

%% Run
source_numbers = [1 6 12];
E1 = zeros(n_designs, n_SNRs, 3);
E2 = zeros(n_designs, n_SNRs, 3);
for kk = 1:3
    doas = linspace(-pi/3, pi/3, source_numbers(kk));
    P = eye(length(doas)) * power_source;
    for dd = 1:n_designs
        design = designs{dd};
        for ii = 1:n_SNRs
            SNR = SNRs(ii);
            power_noise = power_source*10^(-SNR/10);
            e_an = ecov_coarray_music_1d(design, wavelength, doas, power_source,...
                    power_noise, n_snapshots, 'DiagonalsOnly');
            CRB = crb_uc_sto_1d(design, wavelength, doas, P, power_noise, n_snapshots);
            E1(dd, ii, kk) = sum(real(e_an));
            E2(dd, ii, kk) = sum(real(diag(CRB)));
        end
    end
end
efficiency_SNR = E2 ./ E1;

%% Plot
figure;
for kk = 1:3
    subplot(3, 1, kk);
    eff_legends = {};
    markers = {'+', 'o', 'diamond', '*'};
    for dd = 1:n_designs
        hp = plot(SNRs, efficiency_SNR(dd,:,kk), '-', 'Marker', markers{dd}); hold on;
        eff_legends = [eff_legends designs{dd}.name]; 
    end
    hold off;
    xlabel('SNR'); ylabel('\kappa');
    axis([-10 20 0 1]);
    title(sprintf('K = %d', source_numbers(kk)));
    legend(eff_legends{:});
end
set(gcf, 'Position', [200 100 600 800]);
