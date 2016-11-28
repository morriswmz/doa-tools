function pos_perturbed = get_perturbed_positions(pos_nominal, pos_err)
%GET_PERTURBED_POSITIONS Utility function to obtain the perturbed element
%positions.
%Syntax:
%   pos_perturbed = GET_PERTURBED_POSITIONS(pos_nominal, pos_err);
%Inputs:
%   pos_nominal - Nominal element positions.
%   pos_err - Position error matrix.
%Output:
%   pos_perturbed - Perturbed element positions.
a_dim = size(pos_nominal, 1);
p_dim = size(pos_err, 1);
if p_dim < 1 || p_dim > 3
    error('Position error matrix has incorrect number of rows.');
end
pos_perturbed = pos_nominal;
if p_dim < a_dim
    pos_perturbed(1:p_dim,:) = pos_perturbed(1:p_dim,:) + pos_err;
elseif p_dim == a_dim
    pos_perturbed = pos_perturbed + pos_err;
else
    pos_perturbed = pos_perturbed + pos_err(1:a_dim,:);
    pos_perturbed = [pos_perturbed; pos_err(a_dim+1:end,:)];
end
end

