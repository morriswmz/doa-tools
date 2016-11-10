function A = steering_matrix(design, wavelength, doas)
%STEERING_MATRIX Creates the steering matrix for arrays.
% 1D arrays are assumed to be placed along the x-axis.
% 2D arrays are assumed to be placed along the xy-plane.
% For 1D arrays, DOAs correspond to the broadside angle and ranges from
% -pi/2 to pi/2.
% For 2D or 3D arrays, DOAs consist of both azimuth and elevation angles.
%Syntax:
%   A = STEERING_MATRIX(design, wavelength, doas);
%Inputs:
%   design - Array design.
%   wavelength - Wavelength.
%   doas - DOA vector (broadside angles) or matrix (each column represents
%          a pair of azimuth and elevation angles). For 1D arrays, 2D DOAs
%          will be converted to broadside angles. For 2D and 3D arrays, 1D
%          DOAs will be treated as elevation angles.
%Output:
%   A - Steering matrix.
if design.dim == 1
    if ~isvector(doas)
        doas = ae2broad(doas(1,:), doas(2,:));
    else
        doas = reshape(doas, 1, []);
    end
    A = exp(2j*pi/wavelength*(design.element_positions' * sin(doas)));
else
    if isvector(doas)
        doas = reshape(doas, 1, []);
        doas = [zeros(1, length(doas)); doas];
    end
    sin_el = sin(doas(2,:));
    ss = sin_el .* sin(doas(1,:));
    sc = sin_el .* cos(doas(1,:));
    if design.dim == 2
        A = exp(2j*pi/wavelength*(design.element_positions(1,:)' * sc + ...
                design.element_positions(2,:)' * ss));
    elseif design.dim == 3
        A = exp(2j*pi/wavelength*(design.element_positions(1,:)' * sc + ...
                design.element_positions(2,:)' * ss + ...
                design.element_positions(3,:)' * cos(doas(2,:))));
    else
        error('Incorrect array dimension.');
    end
end
end