function [L] = convertM(L1)
    [m,n] = size(L1);
    L = zeros(1,m);
    k = m/n;
    for i=1:n
        for j=1:m
            L((k*(i-1)+1):k*i) = i;
        end
    end
end