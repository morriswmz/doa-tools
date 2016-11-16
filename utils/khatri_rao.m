function Y = khatri_rao(A, B)
%KHATRI_RAO Performs Khatri-Rao product (column-wise Kronecker product).
%Syntax:
%   Y = KHATRI_RAO(A, B)

[r1, c1] = size(A);
[r2, c2] = size(B);
if c1 ~= c2;
    error('Column size mismatch');
end

Y = zeros(r1 * r2, c1);
for i = 1:c1
    Y(:,i) = kron(A(:,i), B(:,i));
end

end

