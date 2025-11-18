
function [normM] = NormMat(M)
% NormMat normalizes each row in a given matrix M.

[m n] = size(M);
sums = sum(M, 2);
sums(find(sums == 0)) = Inf;
factorCol = ones(m, 1) ./ sums;
factor = factorCol(:, ones(1, n));
normM = M.*factor;

return
