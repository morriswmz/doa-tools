function n_source = sn_aic(l, n_sensor, n_snapshot)
%SN_AIC Source number detection using AIC.
%Syntax:
%   n_source = SN_AIC(l, n_sensor, n_snapshot);
%Inputs:
%   l - Eigenvalues of the covariance matrix in descending order.
%   n_sensor - Number of sensors used. Should match the length of l.
%   n_snapshot - Number of snapshots.
%Outputs:
%   n_source - Detected number of sources.
%Remark:
%   AIC is inconsistent, and tends to asymptotically overestimate the
%   number of sources. However, it tends to give a higher probability of a
%   correct decision.
%Reference:
%   H. L. Van Trees, Optimum array processing. New York: Wiley, 2002.

ld = zeros(n_sensor, 1);
for ii = 1:n_sensor
    ld(ii) = sn_ld(l, n_sensor, ii - 1, n_snapshot) + ...
        (ii - 1)*(2*n_sensor - ii + 1);
end
[~, n_source] = min(ld);
n_source = n_source - 1;

end

