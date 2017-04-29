function [blk, At, b] = sdpt3_create_hpsd_blk(n)
%SDPT3_CREATE_HPSD_BLK Creates a Hermitian PSD matrix variable for the
%SDPT3 solver.
%Note that
%   W \succeq 0 \iff U = [ Re(W)  Im(W) ] \succeq 0
%                      [ -Im(W) Re(W) ]
%We declare a 2nx2n symmetric matrix variable U and enforce additional
%constraints according to the structure of U:
%   1. U(i+n,j+n) = U(i,j) \forall 1 <= i <= j <= n
%   2. U(i,j+n) = -U(j,i+n) \forall 1 <= i <= j <= n
%Syntax:
%   [blk, At, b] = SDPT3_CREATE_HPSD_BLK(n);
%Inputs:
%   n - Dimension of the matrix.
%Outputs:
%   blk - Block definition: {'s', 2*n}.
%   At - Transposed constraint matrix (as a submatrix of the full
%        constraint matrix).
%   b - Constraint value vector (as a subvector of the full constraint
%       value vector).
blk = {'s', 2*n};
n_c = n * (n + 1);
n_var = n * (2 * n + 1);
n_sub = n * (n + 1) / 2;
n_nz = 2 * n_c - n;
% triplets for the sparse constraint matrix
idx_row = zeros(n_nz, 1);
idx_col = zeros(n_nz, 1);
vals = zeros(n_nz, 1);
i_c = 0;
i_nz = 0;
% constraints are expressed by
%   A svec(W) = b
%   W(i,j) --> svec(W)(i + j(j-1)/2)
% Re(W) = Re(W)
for jj = 1:n
    for ii = 1:jj
        i_c = i_c + 1;
        i_nz = i_nz + 1;
        idx_row(i_nz) = i_c;
        idx_col(i_nz) = ii + jj*(jj-1)/2;
        vals(i_nz) = 1;
        i_nz = i_nz + 1;
        idx_row(i_nz) = i_c;
        idx_col(i_nz) = ii + n + (jj + n)*(jj + n - 1)/2;
        vals(i_nz) = -1;
    end
end
% Im(W) = -Im(W)^T
% diagonals should be zeros
for ii = 1:n
    i_c = i_c + 1;
    i_nz = i_nz + 1;
    idx_row(i_nz) = i_c; 
    idx_col(i_nz) = ii + (ii + 2*n)*(ii - 1)/2 + n_sub;
    vals(i_nz) = 1;
end
% skew symmetry
for jj = 2:n
    for ii = 1:jj-1
        i_c = i_c + 1;
        i_nz = i_nz + 1;
        idx_row(i_nz) = i_c;
        idx_col(i_nz) = ii + (jj + 2*n)*(jj - 1)/2 + n_sub;
        vals(i_nz) = 1;
        i_nz = i_nz + 1;
        idx_row(i_nz) = i_c;
        idx_col(i_nz) = jj + (ii + 2*n)*(ii - 1)/2 + n_sub;
        vals(i_nz) = 1;
    end
end
At = sparse(idx_row, idx_col, vals, n_c, n_var)';
b = sparse(n_c, 1);
end

