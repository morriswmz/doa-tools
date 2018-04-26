% This script demonstrate that even if several sparse linear arrays
% share the same virtual ULA, they exhibit different performances under
% different SNR settings and different numbers of sources.
% This script will produce results similar to Fig. 1 in the following
% paper:
% * M. Wang, Z. Zhang, and A. Nehorai, "Performance analysis of
%   coarray-based MUSIC and the Cramér-Rao bound," in 2017 IEEE
%   International Conference on Acoustics, Speech and Signal Processing
%   (ICASSP), 2017, pp. 3061-3065.


clear(); close all;

wavelength = 1; % normalized wavelength
d_0 = wavelength / 2;
designs = {...
    design_array_1d('nested', [5  6], d_0, 'Nested (5, 6)') ...
    design_array_1d('nested', [2 12], d_0, 'Nested (2, 12)') ...
    design_array_1d('nested', [3  9], d_0, 'Nested (3, 9)') ...
    design_array_1d('nested', [1 18], d_0, 'Nested (1, 18)') ...
};
n_designs = length(designs);

power_source = 1;
n_snaphots = 1000;

n_grid = 20;
SNRs = linspace(-10, 20, n_grid);

doas1 = linspace(-pi/3, pi/3, 8);

MSEs_SNR_ana1 = zeros(n_designs, n_grid);
for dd = 1:n_designs
    design = designs{dd};
    A = steering_matrix(design, wavelength, doas1);
    for ii = 1:n_grid
        power_noise = power_source*10^(-SNRs(ii)/10);
        MSEs_SNR_ana1(dd, ii) = mean(ecov_coarray_music_1d(design, wavelength, ...
                doas1, power_source, power_noise, n_snaphots, 'DiagonalsOnly'));
    end
end

doas2 = linspace(-pi/3, pi/3, 20);

MSEs_SNR_ana2 = zeros(n_designs, n_grid);
for dd = 1:n_designs
    design = designs{dd};
    A = steering_matrix(design, wavelength, doas2);
    for ii = 1:n_grid
        power_noise = power_source*10^(-SNRs(ii)/10);
        MSEs_SNR_ana2(dd, ii) = mean(ecov_coarray_music_1d(design, wavelength, ...
                doas2, power_source, power_noise, n_snaphots, 'DiagonalsOnly'));
    end
end

figure;
subplot(1,2,1);
semilogy(SNRs, rad2deg(sqrt(MSEs_SNR_ana1)));
xlabel('SNR (dB)'); ylabel('RMSE (deg)'); grid on;
legend(arrayfun(@(x) x{1}.name, designs, 'UniformOutput', false));
title(sprintf('K = %d', length(doas1)));
subplot(1,2,2);
semilogy(SNRs, rad2deg(sqrt(MSEs_SNR_ana2)));
xlabel('SNR (dB)'); ylabel('RMSE (deg)'); grid on;
legend(arrayfun(@(x) x{1}.name, designs, 'UniformOutput', false));
title(sprintf('K = %d', length(doas2)));
