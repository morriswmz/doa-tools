function [est, resolved] = find_doa_est_1d(dg, y, n)
%FIND_DOA_EST_1D Identifies DOA estimates by finding largest peaks in the
%give spectrum.
%Inputs:
%   dg - DOA grid.
%   y - Computed spectrum. Should have the same length as dg.
%   n - Expected number of DOAs.
%Output:
%   est - A row vector of estimated DOAs. If the number of peaks is
%         fewer than the number of DOAs expected, a empty vector is
%         returned.
%   resolved - A boolean indicating whether the direction finding is
%              successful.
[peaks, idx_est] = findpeaks(y);
if length(idx_est) == n
    % juse fine
    est = dg(idx_est);
    resolved = true;
elseif length(idx_est) > n
    % more than expected, choose largest ones
    [~, sort_idx] = sort(peaks);
    est = dg(idx_est(sort_idx(end-n+1:end)));
    est = sort(est);
    resolved = true;
else
    % failed to detect
    est = [];
    resolved = false;
end
end

