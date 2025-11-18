
function [logStr,newM,Ivals,SortedFeatureIndices] = FS (logStr,M,newY,fsMod)
% FS applies Feature Selection using Information Gain to the COLUMNS of M, where:
%
%   - M is a count matrix in which M(i, j) is, e.g., the number of
%     occurrences of the word x_i in the document y_j
%   - newY is the number of requested features (after the FS)
%   - fsMod defines the FS mod. Choose 1 to sort by I(X;y), and 2 to sort
%     by I(X;y)/p(y) (i.e., ignoring feature prior)
%
% FS returns the following output:
%
%   - newM is the input count matrix with only the best newY feature columns
%   - Ivals gives the contribution of every y to I(X;Y) (divided by p(y) if fsMod==2)
%   - SortedFeatureIndices gives the sorted Y indices

Pxy = M./sum(sum(M));
Px_y = NormMat(M')';
Px = sum(Pxy,2);
Py = sum(Pxy,1);
[X Y] = size(Pxy);
I = zeros(1,Y);

fprintf ('Estimating information gain for every y in Y...\n');
for y=1:Y
    if mod(y,1000)==0
        fprintf ('Done %d (out of %d)...\n',y,Y);
    end
    inds = find(Px_y(:,y));
    if fsMod==1 % sort by I(X;y)
        I(y) = Py(y)*sum(Px_y(inds,y).*MyLog(Px_y(inds,y)./Px(inds)));
    else        % sort by I(X;y)/p(y)
        I(y) = sum(Px_y(inds,y).*MyLog(Px_y(inds,y)./Px(inds)));
    end
end

Ivals = I;
[tmp SortedFeatureIndices] = sort(-Ivals);

% Build new matrix using the best newY features (columns) only
newM = M(:,SortedFeatureIndices(1:newY));

% Check for empty rows
Sum_x = sum(newM,2);
zero_rows_inds = find(Sum_x==0);
if any(zero_rows_inds)
    warning ('After FS: %d rows of newM contain no values\n',length(zero_rows_inds));
end

return
