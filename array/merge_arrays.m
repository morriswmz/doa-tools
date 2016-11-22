function design = merge_arrays(varargin)
%MERGE_ARRAYS Merges multiple arrays into one array.
%Syntax:
%   design = merge_arrays(design1, design2, ...);
if nargin == 1
    design = varargin{1};
elseif nargin == 2
    % merge two
    d1 = varargin{1};
    d2 = varargin{2};
    new_dim = max(d1.dim, d2.dim);
    can_have_idx = isfield(d1, 'element_spacing') && isfield(d2, 'element_spacing') && ...
        d1.element_spacing == d2.element_spacing && ...
        isfield(d1, 'element_indices') && isfield(d2, 'element_indices');
    if d1.dim < new_dim
        d1.element_positions = [d1.element_positions; zeros(new_dim - d1.dim, d1.element_count)];
    end
    if d2.dim < new_dim
        d2.element_positions = [d2.element_positions; zeros(new_dim - d2.dim, d2.element_count)];
    end
    new_positions = [d1.element_positions d2.element_positions];
    new_positions = unique(new_positions', 'rows')';
    design.type = 'custom';
    design.dim = new_dim;
    design.name = 'Merged array';
    design.element_positions = new_positions;
    design.element_count = size(new_positions, 2);
    if can_have_idx
        design.element_indices = new_positions / d1.element_spacing;
        design.element_spacing = d1.element_spacing;
    end
elseif nargin == 3
    design = merge_arrays(varargin{1}, merge_arrays(varargin{2}, varargin{3}));
elseif nargin == 4
    design = merge_arrays(merge_arrays(varargin{1}, varargin{2}), merge_arrays(varargin{3}, varargin{4}));
else
    mid = floor(nargin / 2);
    design = merge_arrays(merge_arrays(varargin{1:mid}), merge_arrays(varargin{mid+1:end}));
end
end

