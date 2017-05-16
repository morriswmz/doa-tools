function design = design_array_1d(type, varargin)
%DESIGN_ARRAY_1D Designs a 1D array (linear array).
%Syntax:
%   design = DESIGN_ARRAY_1D('type', ...);
%   design = DESIGN_ARRAY_1D('custom', element_positions[, element_spacing, name]);
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
        design = nested_1d(varargin{:});
    case 'mra'
        design = mra_1d(varargin{:});
    case 'custom'
        pos = varargin{1};
        if ~isvector(pos)
            error('Vector expected for element positions.');
        end
        design = struct();
        design.type = 'custom';
        design.dim = 1;
        design.element_positions = pos;
        design.element_count = length(pos);
        if nargin <= 2
            warning('Element spacing is not specified. Element indices will not be available.');
        else
            design.element_spacing = varargin{2};
            design.element_indices = pos / varargin{2};
        end
        if nargin <= 3
            design.name = 'Custom array';
        else
            design.name = varargin{3};
        end
    otherwise
        error('Unknown array type ''%s''.', type);
end
end

