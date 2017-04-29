function J = svec2vec_conversion_matrix(n)
%SVEC2VEC_CONVERSION_MATRIX Creates a conversion matrix J such that for any
%n dimensional symmetric matrix X, vec(X) = J svec(X).
%Syntax:
%   J = SVEC2VEC_CONVERSION_MATRIX(n);
%Input:
%   n - Dimension.
%Output:
%   J - Conversion matrix of dimension n^2 x n(n+1)/2.
triplets = zeros(n*n, 3);
for jj = 1:n
    for ii = 1:jj-1
        idx_vec = ii + n*(jj - 1);
        triplets(idx_vec, 1) = idx_vec;
        triplets(idx_vec, 2) = ii + jj*(jj - 1)/2;
        triplets(idx_vec, 3) = 1/sqrt(2);
    end
    idx_vec = jj + n*(jj - 1);
    triplets(idx_vec, 1) = idx_vec;
    triplets(idx_vec, 2) = jj + jj*(jj - 1)/2;
    triplets(idx_vec, 3) = 1;
    for ii = jj+1:n
        idx_vec = ii + n*(jj - 1);
        triplets(idx_vec, 1) = idx_vec;
        triplets(idx_vec, 2) = jj + ii*(ii - 1)/2;
        triplets(idx_vec, 3) = 1/sqrt(2);
    end
end
J = spconvert(triplets);
end

