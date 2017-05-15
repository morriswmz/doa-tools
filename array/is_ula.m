function flag = is_ula(design)
%IS_ULA Checks if the given array design is ULA or not.
%Syntax:
%   flag = IS_ULA(design);
%   flag = IS_ULA([0 1 2 3]);
%Input:
%   design - An array design struct or a vector indicating element
%            positions.
%Output:
%   flag - Whether the given array design is ULA.
tolerance = 1e-10;
if isstruct(design)
    if design.dim > 1
        flag = false;
    else
        if design.element_count == 1
            flag = true;
        else
            diffs = diff(design.element_positions);
            flag = all(diff(diffs) / diffs(1) < tolerance);
        end
    end
else
    if ~isvector(design)
        flag = false;
    else
        n = length(design);
        if n == 1
            flag = true;
        else
            diffs = diff(design);
            flag = all(diff(diffs) / diffs(1) < tolerance);
        end
    end
end
end

