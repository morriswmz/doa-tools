function [Rv, dv] = virtual_ula_cov_1d(design, R, mode)
%VIRTUAL_ULA_COV_1D Constructs the augmented covariance matrix of the
% central ULA part of the difference coarray, using either spatial
% smoothing or direct augmentation.
%Syntax:
%Inputs:
%   design - Array design.
%   R - Original sample covariance matrix.
%   mode - Augmentation mode, either 'ss' or 'da'. 'ss' mode uses
%          spatial smoothing and ensures positive-semidefiniteness.
%          'da' mode uses direct augmentation which does not ensure
%          positive-semidefiniteness.
%Outputs:
%   Rv - Augmented covariance matrix.
%   dv - Virtual array design corresponding to Rv.
if nargin <= 2
    mode = 'ss';
end
% extract full difference coarray
[diffs, ~, index_map] = weight_function_1d(design);
% find largest central ULA
w = 0;
zero_idx = find(diffs == 0);
n_diff = length(diffs);
while zero_idx + w < n_diff && zero_idx - w > 1
    w = w + 1;
    if (diffs(zero_idx+w) - diffs(zero_idx+w-1) ~= 1) || (diffs(zero_idx-w+1) - diffs(zero_idx-w)) ~= 1
        w = w - 1;
        break;
    end
end
z = zeros(2*w+1,1);
jj = 1;
for ii=zero_idx-w:zero_idx+w
    z(jj) = mean(R(index_map{ii}));
    jj = jj + 1;
end
Rv = zeros(w + 1);
switch lower(mode)
    case 'ss'
        for ii = 1:w + 1
            Rv = Rv + z(ii:ii + w) * z(ii:ii + w)';
        end
        Rv = Rv / (w + 1);
    case 'da'
        for ii = 1:w + 1
            Rv(:,end-ii+1) = z(ii:ii+w);
        end
    otherwise
        error('Unknown mode "%s".', mode);
end
dv = ula_1d(w + 1, design.element_spacing, 'Virtual ULA');
end

