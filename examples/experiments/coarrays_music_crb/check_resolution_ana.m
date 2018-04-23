function flag = check_resolution_ana(design, wavelength, doa1, doa2, power_source, noise_var, n_snapshots, lambda)
%CHECK_RESOLUTION_ANA Checks if two sources are resolvable analytically with
%SS-MUSIC (or DA-MUSIC) using the results in the following paper:
% * M. Wang and A. Nehorai, "Coarrays, MUSIC, and the Cramer-Rao Bound,"
%   IEEE Trans. Signal Process., vol. 65, no. 4, pp. 933-946, Feb. 2017.

if nargin == 7
    lambda = 1;
end

% ensure doa1 < doa2
if doa1 > doa2
    tmp = doa1;
    doa1 = doa2;
    doa2 = tmp;
end

% compute P
if isscalar(power_source)
    P = power_source*eye(2);
else
    if length(power_source) ~= 2
        error('Length of source powers does not match the number of sources.');
    end
    P = diag(power_source);
end

se = ecov_coarray_music_1d(design, wavelength, [doa1 doa2], P, noise_var, n_snapshots, 'DiagonalsOnly');
% resolution criterion
flag = doa2 - doa1 > lambda*(sqrt(se(1)) + sqrt(se(2)));

end

