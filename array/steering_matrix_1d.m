function A = steering_matrix_1d(design, wavelength, doas)
%STEERING_MATRIX_1D Creates the steering matrix for 1D arrays.
%Syntax:
%   A = steering_matrix_1d(design, wavelength, doas);
%Inputs:
%   design - Array design.
%   wavelength - Wavelength.
%   doas - DOA vector.
%Output:
%   A - Steering matrix.
A = exp(-2j*pi/wavelength*(design.element_positions * sin(reshape(doas, 1, []))));
end

