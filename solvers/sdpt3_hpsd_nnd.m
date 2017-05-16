function [R, info] = sdpt3_hpsd_nnd(R_est, lambda, varargin)
%SDPT3_HPSD_NND Nuclear norm denoising of Hermitian PSD matrix. By default,
%this solves the following optimization problem.
%   (P1) min_R || R ||_*
%        s.t.  || R - \hat{R} ||_F <= \lambda
%              R \succeq 0, R^H = R
%The penalized form is given by
%   (P2) min_R || R - \hat{R} ||_F^2 + \lambda || R ||_*
%        s.t.  R \succeq 0, R^H = R
%Syntax:
%   [R, info] = SDPT3_HPSD_NND(R_est, lambda, ...);
%Inputs:
%   R_est - Noisy estimate of the Hermitian PSD matrix R.
%   lambda - Regularization parameter.
%   ... - Options:
%       'Formulation' - Either 'Constrained' (P1) or 'Penalized' (P2).
%                       Default value is 'Constrained'.
%       'Verbose' - Set to true to enable detailed output from the SDPT3
%                   solver. Default value is false.
%       'SQLPOptions' - If present, will override the default options for
%                       the SDPT3 solver. This also overrides the print
%                       level set by the 'Verbose' option.
%       'W' - n^2 x n^2 weighting matrix for the norm computation. Can be
%                       complex but must be PSD. If provided, the Frobenius
%                       norm term will be replaced with
%                           || vec(R) - vec(\hat{R}) ||_W
%Outputs:
%   R - Denoised version of R_est.
%   info - Solution information provided by SDPT3. See sqlp.m for details.

use_penalized_formulation = false;
verbose = false;
sqlp_options = [];
W = [];
for ii = 1:2:length(varargin)
    option_name = varargin{ii};
    option_value = varargin{ii+1};
    switch lower(option_name)
        case 'formulation'
            switch lower(option_value)
                case 'constrained'
                    use_penalized_formulation = false;
                case 'penalized'
                    use_penalized_formulation = true;
                otherwise
                    error('Unknown formulation ''%s''.', option_value);
            end
        case 'verbose'
            verbose = option_value;
        case 'sqlpoptions'
            sqlp_options = option_value;
        case 'w'
            W = option_value;
        otherwise
            error('Unknown option ''%s''.', option_name);
    end
end

n = size(R_est, 1);
n_blk = 2;
blk = cell(n_blk,2);
At = cell(n_blk,1);
C = cell(n_blk,1);

% Because R is Hermitian, || R ||_* = tr(R). For a general matrix R,
% nuclear norm minimization problem can be casted into the following
% problem:
%   min tr(W_1) + tr(W_2)
%   s.t. U = [ W_1 R   ]
%            [ R   W_2 ] \succeq 0
%        W_1 \succeq 0, W_1 = W_1^H
%        W_2 \succeq 0, W_2 = W_2^H
%        R \succeq 0, R = R^H
% which is much complex to implement for SDPT3.
% Here we only consider R to be Hermitian PSD. To tackle the complex
% matrices, we use the following property:
%   R \succeq 0 <=> [ Re(R) -Im(R) ] \succeq 0
%                   [ Im(R)  Re(R) ] 
% Also, under svec indexing,
%   X(i,j) -> svec(X)(i + j(j-1)/2) (ignoring the scaling factor)

% R: 2n x 2n
[blk_R, At_R, ~] = sdpt3_create_hpsd_blk(n);
blk(1,:) = blk_R;
n_var_s = n*(2*n + 1);
n_c_s = size(At_R, 2);
% for the penalized formulation, we convert the optimization problem into
%   min_R t + \lambda || R ||_*
%   s.t. || R - \hat{R} ||_F <= z
%        [ 1 z ] \succeq 0
%        [ z t ]
%        R \succeq 0, R^H = R
% we need to add an additional 2x2 block
if use_penalized_formulation
    blk{1,2} = [blk{1,2} 2];
    At_R = [At_R; sparse(3, n_c_s)];
    n_var_s = n_var_s + 3;
end

% SOC constraint
% z = [lambda; vec(real(R)); vec(imag(R)) excluding diagonals];
% because B can be non Hermitian
% unfortunately we need additional constraints here
% TODO: optimize when B is Hermitian using the following results:
% || X - B ||_F^2 = || X - (B + B^H)/2 ||_F^2 +
%   0.5 (|| B ||_F^2 - tr(real(BB)))
n_var_q = 2*n*n + 1;
blk{2,1} = 'q'; blk{2,2} = n_var_q;

% convert W
% W: 2n^2 x 2n^2
if isempty(W)
    W = speye(2*n*n);
else
    W_re = real(W);
    W_im = imag(W);
    W = [W_re -W_im;W_im W_re];
    W = chol(W);
end
% convert R_est
r_est = [real(R_est(:)); imag(R_est(:))];
% now the constraint becomes
% - W [vec(real(R)); vec(imag(R))] + z(2:end) = - W r_est
r_est = W * r_est;
% generate mapping from real representation of
%   svec ( [ Re(R) -Im(R) ] )
%        ( [ Im(R)  Re(R) ] )
% to
%   [vec(real(R)); vec(imag(R))]
J_re = [svec2vec_conversion_matrix(n) sparse(n*n, n_var_s - n*(n + 1)/2)];
triplets = zeros(n*n, 3);
for jj = 1:n
    for ii = 1:jj-1
        idx_vec = ii + n*(jj - 1);
        triplets(idx_vec, 1) = idx_vec;
        triplets(idx_vec, 2) = ii + (jj + n)*(jj + n - 1)/2;
        triplets(idx_vec, 3) = -1/sqrt(2);
    end
    idx_vec = jj + n*(jj - 1);
    triplets(idx_vec, 1) = idx_vec;
    triplets(idx_vec, 2) = jj + (jj + n)*(jj + n - 1)/2;
    triplets(idx_vec, 3) = -1/sqrt(2);
    for ii = jj+1:n
        idx_vec = ii + n*(jj - 1);
        triplets(idx_vec, 1) = idx_vec;
        triplets(idx_vec, 2) = jj + (ii + n)*(ii + n - 1)/2;
        triplets(idx_vec, 3) = 1/sqrt(2);
    end
end
J_im = sparse(triplets(:,1), triplets(:,2), triplets(:,3), n*n, n_var_s);
At_Rz = -W*[J_re; J_im];
At_Rz = sparse(At_Rz');

% construct full constraint matrices
if use_penalized_formulation
    % note that we need to deal with the scaling factor sqrt(2) resulting
    % from svec here
    At{1} = [At_R ...
        sparse(n_var_s - 1, 1, -1/sqrt(2), n_var_s, 1) ...
        At_Rz ...
        sparse(n_var_s - 2, 1, 1, n_var_s, 1)];
    At{2} = [sparse(n_var_q, n_c_s) speye(n_var_q) sparse(n_var_q, 1)];
else
    At{1} = [At_R sparse(n_var_s, 1) At_Rz];
    At{2} = [sparse(n_var_q, n_c_s) speye(n_var_q)];
end
% assemble b, note the minus sign before r_est
if use_penalized_formulation
    b = [sparse(n_c_s,1); 0; -r_est; 1];
else
    b = [sparse(n_c_s,1); lambda; -r_est];
end

% objective function
C{2} = sparse(n_var_q, 1);
if use_penalized_formulation
    % \lambda * tr(R) + t
    C1 = 0.5 * lambda * speye(2*n + 2, 2*n + 2);
    C1(2*n+1, 2*n+1) = 0;
    C1(2*n+2, 2*n+2) = 1;
    C{1} = C1;
else
    C{1} = 0.5 * speye(2*n, 2*n); % tr(R)
end

% call SDPT3
if isempty(sqlp_options)
    sqlp_options = sqlparameters;
    if ~verbose
        sqlp_options.printlevel = 0;
    end
end
[~,solution,~,~,info,~] = sqlp(blk, At, C, b, sqlp_options);
R = full(solution{1});
R = R(1:n,1:n) - 1j*R(1:n,n+1:2*n);
R = 0.5*(R + R');

end

