function v = logdeth(A, force_hermitian)
%LOGDETH Log-determinant function for Hermitian symmetric matrices. This
%function uses Cholesky decomposition. If the decomposition fails, inf will be
%returned. This input matrix is assumed Hermitian.
if nargin == 1 || force_hermitian
    A = 0.5*(A+A');
end
[R, p] = chol(A);
if p > 0
    v = inf;
else
    v = sum(log(diag(R))) * 2;
end
end