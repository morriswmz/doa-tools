function s = extend_struct(s, varargin)
%EXTEND_STRUCT Extends a structure (matrix) with new fields.
%Syntax:
%   s = EXTEND_STRUCT(s, ['Prop1', val_1, ...]);
%   s = EXTEND_STRUCT(s, s_1, s_2, ...);
%Inputs:
%   s - The original structure (matrix). Can be empty.
%   s_1, s_2, ... - Structures whose fields will be copied to s. Each of
%                   them must be either a scalar, or a matrix sharing the
%                   same dimensions as s. If s_i is a scalar, then
%                       s(k).* = s_i.* for all k
%                   If s_i is a matrix, then
%                       s(k).* = s_i(k).* for all k
%   'Prop1', val_1, ... - Name-value pairs describing the new
%                   fields. By default, the following assigment is
%                   performed:
%                       s(k).Propi = val_i for all k
%                   If val_i is a matrix, and you want the following
%                   behavior:
%                       s(k).Propi = val_i(k) for all k
%                   append '[]' to the field name such that 'Propi'
%                   becomes 'Propi[]'.
%                   The same applies when val_i is a cell array. In this
%                   case, append '{}' instead.
%Output:
%   s - Extended structure (matrix).
%Note:
%   Existing fields can be overriden, whose values will be determined by
%   the last assignment. For example:
%       EXTEND_STRUCT([], 'v', 1, 'v', 2)
%   will return a structure s with s.v = 2;
%Examples:
%   % add a new field
%   s = struct('value', 12); % {"value": 12}
%   s1 = EXTEND_STRUCT(s, 'correct', true); % {"value": 12, "correct": true}
%   % overrid an existing field
%   s2 = struct('value', 9); % {"value": 9}
%   s3 = EXTEND_STRUCT(s, s2); % {"value": 9}
%   % use the '[]' affix
%   s4 = EXTEND_STRUCT(repmat(struct(), 1, 2), 'value[]', [3 4]);
%   % [{"value": 3}, {"value": 4}]
if nargin <= 1
    s = struct();
    return;
end
if isempty(s)
    s = struct();
end
if ischar(varargin{1})
    % name-value pairs
    if mod(nargin - 1, 2) ~= 0
        error('Incomplete name-value pairs.');
    end
    if isscalar(s)
        for ii = 1:2:nargin - 1
            field_name = varargin{ii};
            if has_matrix_affix(field_name) || has_cell_affix(field_name)
                error('Using ''[]'' or ''{}'' affix is not allowed when the input structure is a scalar.');
            end
            s.(field_name) = varargin{ii+1};
        end
    else
        for ii = 1:2:nargin - 1
            field_name = varargin{ii};
            if has_matrix_affix(field_name)
                field_name = field_name(1:end - 2);
                vals = varargin{ii + 1};
                if any(size(vals) ~= size(s))
                    error('The size of the value matrix does not match that of the structure.');
                end
                % Note: in this case, for loop is faster because num2cell
                % conversion is slow.
                n_el = numel(s);
                for nn = 1:n_el
                    s(nn).(field_name) = vals(nn);
                end
            else
                % Note: the following destructing assigment is much faster
                % than for loops. MATLAB magic!
                if has_cell_affix(field_name)
                    field_name = field_name(1:end - 2);
                    vals = varargin{ii + 1};
                    if any(size(vals) ~= size(s))
                        error('The size of the value cell does not match that of the structure.');
                    end
                else
                    vals = repmat(varargin(ii + 1), numel(s), 1);
                end
                [s.(field_name)] = vals{:};
            end
        end
    end
elseif isstruct(varargin{1})
    % a struct
    for ss = 1:nargin - 1
        s2 = varargin{ss};
        new_fields = fieldnames(s2);
        if isscalar(s2)
            if isscalar(s)
                for ii = 1:length(new_fields)
                    s.(new_fields{ii}) = s2.(new_fields{ii});
                end
            else
                for ii = 1:length(new_fields)
                    % see the note above
                    vals = repmat({s2.(new_fields{ii})}, numel(s), 1);
                    [s.(new_fields{ii})] = vals{:};
                end
            end
        else
            if any(size(s2) ~= size(s))
                error('Structure matrix dimensions mismatch.');
            end
            for ii = 1:length(new_fields)
                % see the note above
                vals = {s2.(new_fields{ii})};
                [s.(new_fields{ii})] = vals{:};
            end
        end
        
    end
else
    error('Invalid syntax.');
end
end

function flag = has_matrix_affix(field_name)
flag = length(field_name) > 2 && strcmp(field_name(end - 1:end), '[]');
end

function flag = has_cell_affix(field_name)
flag = length(field_name) > 2 && strcmp(field_name(end - 1:end), '{}');
end

