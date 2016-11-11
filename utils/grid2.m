function g = grid2(xmin, xmax, ymin, ymax, nx, ny)
%GRID2 Creates a 2D parameter grid.
%Syntax:
%   g = GRID2(xmin, xmax, ymin, ymax, nx[, ny]);
%Inputs:
%   xmin - Minimum x value (inclusive).
%   xmax - Maximum x value (inclusive).
%   ymin - Minimum y value (inclusive).
%   ymax - Maximum y value (inclusive).
%   nx - Number of grid point along x.
%   ny - (Optional) number of grid points along y. If omitted, will be set
%        to nx.
%Output:
%   g - A (nx x ny) x 2 parameter grid arranged as
%           x_1 x_1 ... x_1  ... x_nx x_nx ... x_nx
%           y_1 y_2 ... y_ny ... y_1  y_2  ... y_ny
if nargin <= 5 || isempty(ny)
    ny = nx;
end
xg = linspace(xmin, xmax, nx);
yg = linspace(ymin, ymax, ny);
g = [reshape(repmat(xg, ny, 1), 1, []); repmat(yg, 1, nx)]; 
end

