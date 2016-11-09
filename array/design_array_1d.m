function design = design_array_1d(type, varargin)
%DESIGN_ARRAY_1D Designs a 1D array (linear array).
%Syntax:
%   design = DESIGN_ARRAY_1D('type', ...);
%Inputs:
%   type - Type of the array.
%   ... - Specific array design parameters.
%Output:
%   design - Array design.
switch lower(type)
    case {'ula', 'uniform'}
        design = ula_1d(varargin{:});
    case {'coprime', 'co-prime'}
        design = coprime_1d(varargin{:});
    case 'nested'
    case 'mra'
    case 'custom'
    otherwise
        error('Unknown array type "%s".', type);
end
end

