function n_source = sn_mdl(l, n_sensor, n_snapshot)
%SN_MDL Source number detection using MDL.
%Syntax:
%   n_source = SN_MDL(l, n_sensor, n_snapshot);
%Inputs:
%   l - Eigenvalues of the covariance matrix in descending order.
%   n_sensor - Number of sensors used. Should match the length of l.
%   n_snapshot - Number of snapshots.
%Outputs:
%   n_source - Detected number of sources.
%Remark:
%   MDL is consistent.
%Reference:
%   H. L. Van Trees, Optimum array processing. New York: Wiley, 2002.


ld = zeros(n_sensor, 1);
for ii = 1:n_sensor
    ld(ii) = sn_ld(l, n_sensor, ii - 1, n_snapshot) + ...
        0.5*((ii - 1)*(2*n_sensor - ii + 1) + 1)*log(n_snapshot);
end
[~, n_source] = min(ld);
n_source = n_source - 1;

end

