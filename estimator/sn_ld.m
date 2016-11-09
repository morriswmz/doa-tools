function ld = sn_ld(l, n_sensor, n_source, n_snapshot)
%SN_LD Sufficient statistic for source number detection in MDL/AIC.
%   This function is used internally.
%Syntax:
%   ld = SN_LD(l, n_sensor, n_source, n_snapshot);
%Inputs:
%   l - Eigenvalues of the covariance matrix in descending order.
%   n_sensor - Number of sensors used. Should match the length of l.
%   n_source - Number of sources.
%   n_snapshot - Number of snapshots.
%Outputs:
%   ld - Computed value.
%Reference:
%   H. L. Van Trees, Optimum array processing. New York: Wiley, 2002.

diff = n_sensor - n_source;
l = l(n_source+1:end);
ld = n_snapshot * diff * log(sum(l)/diff/(prod(l)^(1/diff)));

end

