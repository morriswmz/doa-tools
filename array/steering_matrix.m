function [A, DA] = steering_matrix(design, wavelength, doas, use_normalized_angles)
%STEERING_MATRIX Creates the steering matrix for arrays.
% 1D arrays are assumed to be placed along the x-axis.
% 2D arrays are assumed to be placed along the xy-plane.
% For 1D arrays, DOAs correspond to the broadside angle and ranges from
% -pi/2 to pi/2. By the above assumptions, the broadside angle is always
% measured from the yz-plane.
% For 2D or 3D arrays, DOAs consist of both azimuth and elevation angles.
%Syntax:
%   A = STEERING_MATRIX(design, wavelength, doas, ...);
%   [A, DA] = STEERING_MATRIX(design, wavelength, doas, ...);
%Inputs:
%   design - Array design. Can optionally contain fields describing the
%            model errors:
%               mutual_coupling : A element_count x element_count matrix
%               gain_errors: Gain error vector of length element_count.
%               phase_errors: Phase error vector of length element_count.
%               position_errors: Position error matrix of size
%                   p_dim x element_count, where p_dim denotes the number
%                   of dimensions of the postion errors. If p_dim is 1,
%                   positioning_errors specifies the perturbation along the
%                   x-axis. If p_dim is 3, positioning_errors specifies the
%                   peturbations along all the three axis.
%   wavelength - Wavelength.
%   doas - DOA vector (broadside angles) or matrix (each column represents
%          a pair of azimuth and elevation angles). For 1D arrays, 2D DOAs
%          will be converted to broadside angles. For 2D and 3D arrays, 1D
%          DOAs will be treated as elevation angles.
%   use_normalized_angles - If set to true, the input DOAs will be treated
%          as normalized angles. Given a angle theta, the normalized angle
%          is given by
%                   2*pi*wavelength*sin(theta) / d_0.
%          In this case, the wavelength parameter is ignored. 
%          Default value is false.
%Output:
%   A - Steering matrix.
%   DA - Derivative matrix, the i-th column of which is the derivative of
%        the i-th column of A with respect to the i-th DOA.
%        Only available to the 1D DOA case.
if nargin < 4
    use_normalized_angles = false;
end
if use_normalized_angles
    if ~isfield(design, 'element_spacing')
        error('The element_spacing field must be specified in the design.');
    end
end
need_derivative = nargout == 2;
if need_derivative && ~isvector(doas)
    error('Derivative matrix is only available for 1D DOA estimation.');
end
% handle position errors
if isfield(design, 'position_errors') && ~isempty(design.position_errors)
    % change design dim here
    design.element_positions = get_perturbed_positions(design.element_positions, design.position_errors);
    design.dim = size(design.element_positions, 1);
end
% compute the steering matrix
DA = [];
if design.dim == 1
    if ~isvector(doas)
        doas = ae2broad(doas(1,:), doas(2,:));
    else
        doas = reshape(doas, 1, []);
    end
    if use_normalized_angles
        pos = design.element_positions' / design.element_spacing;
        A = exp(1j * pos * doas);
        if need_derivative
            DA = bsxfun(@times, (1j * pos), A);
        end
    else
        A = exp(2j*pi/wavelength*(design.element_positions' * sin(doas)));
        if need_derivative
            DA = A .* (2j*pi/wavelength*(design.element_positions' * cos(doas)));
        end
    end
else
    if design.dim ~= 2 && design.dim ~= 3
        error('Incorrect array dimension.');
    end
    if use_normalized_angles
        error('Normalized angles are not supported for 2D/3D arrays yet.');
    end
    if isvector(doas)
        % only broadside angles are provided
        % to be consistent with the 1D case where the array is lied on the
        % x-axis, we assume the sources are on the xy-plane, and the
        % broadside angle here is relative to the y-axis (positive on the
        % right side of the y-axis).
        doas = reshape(doas, 1, []);
        % for 3D arrays, the offsets along z-axis does not matter
        A = exp(2j*pi/wavelength*(design.element_positions(1,:)' * sin(doas) + ...
                design.element_positions(2,:)' * cos(doas)));
        if need_derivative
            DA = A .* (2j*pi/wavelength*(design.element_positions(1,:)' * cos(doas) - ...
                    design.element_positions(2,:)' * sin(doas)));
        end
    else
        % normal azimuth-elevation pairs
        cos_el = cos(doas(2,:));
        cs = cos_el .* sin(doas(1,:));
        cc = cos_el .* cos(doas(1,:));
        if design.dim == 2
            A = exp(2j*pi/wavelength*(design.element_positions(1,:)' * cc + ...
                    design.element_positions(2,:)' * cs));
        elseif design.dim == 3
            A = exp(2j*pi/wavelength*(design.element_positions(1,:)' * cc + ...
                    design.element_positions(2,:)' * cs + ...
                    design.element_positions(3,:)' * sin(doas(2,:))));
        end
    end
end
% handle mutual coupling, gain and phase errors
if isfield(design, 'mutual_coupling') && ~isempty(design.mutual_coupling)
    A = design.mutual_coupling * A;
    if ~isempty(DA)
        DA = design.mutual_coupling * DA;
    end
end
if isfield(design, 'gain_errors') && ~isempty(design.gain_errors)
    A = bsxfun(@times, design.gain_errors(:), A);
    if ~isempty(DA)
        DA = bsxfun(@times, design.gain_errors(:), DA);
    end
end
if isfield(design, 'phase_errors') && ~isempty(design.phase_errors)
    phase_errors = exp(1j*design.phase_errors(:));
    A = bsxfun(@times, phase_errors, A);
    if ~isempty(DA)
        DA = bsxfun(@times, phase_errors, DA);
    end
end
end