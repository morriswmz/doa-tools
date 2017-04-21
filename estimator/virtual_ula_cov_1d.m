function [Rv, dv, z] = virtual_ula_cov_1d(design, R, mode)
%VIRTUAL_ULA_COV_1D Constructs the augmented covariance matrix of the
%positive part of the central ULA part of the difference coarray, using
%either spatial smoothing or direct augmentation.
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
%   z - Virtual measurement vector of the central ULA (includes both the
%       positive part and the negative part).
if nargin <= 2
    mode = 'ss';
end
% extract full difference coarray
[diffs, ~, index_map] = weight_function_1d(design);
% find largest central ULA
[s, m_v, zero_idx] = get_central_ula_size(diffs);
z = zeros(s,1);
jj = 1;
for ii=zero_idx-m_v+1:zero_idx+m_v-1
    z(jj) = mean(R(index_map{ii}));
    jj = jj + 1;
end
Rv = zeros(m_v);
switch lower(mode)
    case 'ss'
        for ii = 1:m_v
            Rv = Rv + z(ii:ii+m_v-1) * z(ii:ii+m_v-1)';
        end
        Rv = Rv / m_v;
    case 'da'
        for ii = 1:m_v
            Rv(:,end-ii+1) = z(ii:ii+m_v-1);
        end
    otherwise
        error('Unknown mode "%s".', mode);
end
dv = ula_1d(m_v, design.element_spacing, 'Virtual ULA');
end
