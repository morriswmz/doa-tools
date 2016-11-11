function [A, dg_rad, dg_display, dg_range] = default_steering_matrix_grid(design, wavelength, n, unit, dim)
%DEFAULT_STEERING_MATRIX_GRID Creates steering matrix for the discretized
%DOA parameter space. Internal use only.
%Inputs:
%   design - Array design or steering matrix generator.
%   wavelength - Wavelength. If design is a function handle, this parameter
%                should be set to [].
%   n - Number of grid points. For 2D DOAs, n can be 2-element vector,
%       specifying the number of grid points for azimuth and elevation
%       angles, respectively.
%   unit - Can be 'radian', 'degree', or 'sin'.
%   dim - (Optional) DOA dimension. If unspecified, will be determined by
%         the array design.
%Outputs:
%   A - Steering matrix.
%   dg_rad - Grid of candidate DOAs converted to radians.
%   dg_display - Grid of candidate DOAs in the specified unit.
%   dg_range - Range of the DOA candidates. For the 1D case, this is an 1x2
%              vector. For the 2D case, this is a 2x2 matrix, where the
%              first row represents the range of the azimuth angle, and the
%              second row represents that of the elevation angle.
if isempty(dim)
    if design.dim > 1
        dim = 2;
    else
        dim = 1;
    end
else
    if dim ~= 1 && dim ~= 2
        error('Incorrect DOA dimension.');
    end
end
if dim == 2
    if isscalar(n)
        n = [n n];
    end
end
switch lower(unit)
    case 'radian'
        if dim == 1
            dg_range = [-pi/2 pi/2];
            dg_rad = linspace(-pi/2, pi/2, n);
            dg_display = dg_rad;
        else
            dg_range = [0 pi;-pi/2 pi/2];
            dg_rad = grid2(0, pi, -pi/2, pi/2, n(1), n(2));
            dg_display = dg_rad;
        end
    case 'degree'
        if dim == 1
            dg_range = [-pi/2 pi/2];
            dg_rad = linspace(-pi/2, pi/2, n);
            dg_display = rad2deg(dg_rad);
        else
            dg_range = [0 pi;-pi/2 pi/2];
            dg_rad = grid2(0, pi, -pi/2, pi/2, n(1), n(2));
            dg_display = rad2deg(dg_rad);
        end
    case 'sin'
        if dim == 1
            dg_range = [-1 1];
            dg_display = linspace(-1, 1, n);
            dg_rad = asin(dg_display);
        else
            dg_range = [-1 1;-1 1];
            dg_display = grid2(-1, 1, -1, 1, n(1), n(2));
            dg_rad = asin(dg_display);
            dg_rad(1,:) = dg_rad(1,:) + pi/2; % azimuth -> [0, pi]
        end
    otherwise
        error('Unknown unit "%s".', unit);
end
% create steering matrix
if ishandle(design) && isempty(wavelength)
    A = design(wavelength, dg_rad);
else
    A = steering_matrix(design, wavelength, dg_rad);
end
end


