function [est, est_idx, resolved] = find_doa_from_spectrum_1d(dg, y, n)
%FIND_DOA_FROM_SPECTRUM_1D Identifies DOA estimates by finding largest
%peaks in the give spectrum.
%Inputs:
%   dg - DOA grid.
%   y - Computed spectrum. Should have the same length as dg.
%   n - Expected number of DOAs.
%Output:
%   est - A row vector of estimated DOAs. If the number of peaks is
%         fewer than the number of DOAs expected, a empty vector is
%         returned.
%   est_idx - A row vector indicating the index of each estimate in dg
%             (i.e. dg(est_idx) = est).
%   resolved - A boolean indicating whether the direction finding is
%              successful.
[peaks, idx_est] = findpeaks(y);
if length(idx_est) == n
    % juse fine
    est = dg(idx_est);
    est_idx = idx_est;
    resolved = true;
elseif length(idx_est) > n
    % more than expected, choose largest ones
    [~, sort_idx] = sort(peaks);
    est_idx = idx_est(sort_idx(end-n+1:end)); 
    est = dg(est_idx);
    [est, sort_idx] = sort(est);
    est_idx = est_idx(sort_idx);
    resolved = true;
else
    % failed to detect
    est = [];
    est_idx = [];
    resolved = false;
end
end

