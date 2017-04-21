function est = refine_grid_estimates(f_obj, grid, est_idx)
%REFINE_GRID_ESTIMATES Refines DOA estimates obtained from a grid.
%Inputs:
%   f_obj - Objective function. Its local minimums identify the DOAs.
%   grid - Grid used for the original estimation.
%   est_idx - Indices (corresponding to the grid) of the original
%             estimates.
%Output:
%   est - Refined estimates.
est = zeros(size(est_idx));
if size(grid, 1) == 1
    % 1d
    n_iter = 10;
    subgrid_size = 10;
    for kk = 1:length(est_idx)
        % k-th DOA
        % init bounds
        if est_idx(kk) > 1
            lb = grid(est_idx(kk) - 1);
        else
            lb = grid(1);
        end
        if est_idx(kk) < length(grid)
            ub = grid(est_idx(kk) + 1);
        else
            ub = grid(end);
        end
        % refine
        for ii = 1:n_iter
            % find minimum over the subgrid
            subgrid = linspace(lb, ub, subgrid_size);
            obj_vals = f_obj(subgrid);
            [~, min_idx] = min(obj_vals);
            % update bounds
            if min_idx > 1
                lb = subgrid(min_idx - 1);
            else
                lb = subgrid(1);
            end
            if min_idx < subgrid_size
                ub = subgrid(min_idx + 1);
            else
                ub = subgrid(end);
            end
        end
        est(kk) = subgrid(min_idx);
    end
else
    % 2d
    error('Not implemented.');
end
end

