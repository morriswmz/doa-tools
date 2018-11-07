function k = sn_sorte(l)
%SN_SORTE Source number detection using SORTE.
%Input:
%   l - Eigenvalues sorted in descending order.
%Output:
%   k - Estimated source number.
%Reference:
%   [1] Z. He, A. Cichocke, S. Xie, and K. Choi, "Detecting the number of
%       clusters in n-way probabilistic clustering," IEEE Trans. Pattern
%       Anal. Mach. Intell., vol. 32, pp. 2006-2021, Nov. 2010.
n = length(l);
if n < 3
    error('At least three eigenvalues required.');
end
dl = l(1:end-1) - l(2:end);
sum_dl = cumsum(dl, 'reverse');
vars = zeros(n - 1, 1);
for kk = 1:n - 1
    vars(kk) = sum((dl(kk:end) - sum_dl(kk)/(n - kk)).^2)/(n - kk);
end
s = zeros(n - 2, 1);
for kk = 1:n - 2
    if abs(vars(kk)) < eps
        s(kk) = inf;
    else
        s(kk) = vars(kk + 1)/vars(kk);
    end
end
[~, k] = min(s);
end

