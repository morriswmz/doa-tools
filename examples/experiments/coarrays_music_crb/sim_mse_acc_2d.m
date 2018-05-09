% Checks the accuracy of the proposed MSE express w.r.t. different SNRs and
% number of snapshots.
% This script will generate Fig. 2 in the following paper:
% * M. Wang and A. Nehorai, "Coarrays, MUSIC, and the Cramer-Rao Bound,"
%   IEEE Trans. Signal Process., vol. 65, no. 4, pp. 933-946, Feb. 2017.


clear(); close all;

%% Setup
%  =====

use_simple_rss = true;
wavelength = 1; % normalized
d_0 = wavelength / 2;
doas = linspace(-pi*6/16, pi*5/16, 11);
n_doas = length(doas);
power_source = 1;

default_SNR = 0;
default_n_snapshots = 500;
default_power_source = 1;

n_param_n_snapshots = 20;
n_param_SNR = 20;

min_n_snapshots = 10;
max_n_snapshots = 1000;
min_SNR = -20;
max_SNR = 20;

n_repeats = 5000;

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

% params
params_n_snapshots = round(linspace(min_n_snapshots, max_n_snapshots, n_param_n_snapshots));
params_SNR = linspace(min_SNR, max_SNR, n_param_SNR);
rv_modes = {'DA', 'SS'};

%% Run
% ==========
results_em = zeros(n_designs, 2, n_param_n_snapshots, n_param_SNR);
results_an = zeros(n_designs, 2, n_param_n_snapshots, n_param_SNR);

% You can adjust the parfor statements here. Just remember that
% progressbar cannot be called within parfor.
progressbar('reset', n_designs * 2 * n_param_n_snapshots);
for dd = 1:n_designs
    design = designs{dd};
    for ss = 1:2
        rv_mode = rv_modes{ss};
        for ii = 1:n_param_n_snapshots
            n_snapshots = params_n_snapshots(ii);
            parfor jj = 1:n_param_SNR
                SNR = params_SNR(jj);
                power_noise = 10^(-SNR/10) * min(power_source);
                % run simulations
                cur_mse = 0;
                for rr = 1:n_repeats                
                    [~, R] = snapshot_gen_sto(design, doas, wavelength, n_snapshots, power_noise, power_source);
                    Rv = virtual_ula_cov_1d(design, R, rv_mode);
                    sp = rmusic_1d(Rv, n_doas, 2*pi*d_0/wavelength);
                    cur_mse = cur_mse + sum((sp.x_est - doas).^2) / n_doas;
                end
                results_em(dd, ss, ii, jj) = cur_mse / n_repeats;
                e_an = ecov_coarray_music_1d(design, wavelength, doas, power_source, power_noise, n_snapshots, 'DiagonalsOnly');
                results_an(dd, ss, ii, jj) = mean(e_an);
            end
            progressbar('advance');
        end
    end
end
progressbar('end');

%% Plot
%  ====
rel_errors = abs(results_an - results_em) ./ results_em;
figure;
for dd = 1:n_designs
    nn = 2*(dd-1)+1;
    subplot(n_designs, 2, nn);
    imagesc(params_SNR, params_n_snapshots, squeeze(rel_errors(dd,1,:,:)));
    set(gca,'YDir','normal')
    colorbar;
    colormap gray;
    colormap(flipud(colormap));
    if dd == n_designs
        xlabel('SNR (dB)');
    end
    ylabel('Number of snapshots');
    title([designs{dd}.name ' (DA-MUSIC)']);
    hold on;
    
    nn = 2*dd;
    subplot(n_designs, 2, nn);
    imagesc(params_SNR, params_n_snapshots, squeeze(rel_errors(dd,2,:,:)));
    set(gca,'YDir','normal')
    colorbar;
    colormap gray;
    colormap(flipud(colormap));
    set(gca, 'Clim', [0 1]);
    if dd == n_designs
        xlabel('SNR (dB)');
    end
    title([designs{dd}.name ' (SS-MUSIC)']);
end
hold off;
