function [co,new_t] =MinSumCosts(cellCosts,prm,Weights)
% sums=cellCosts{1};
% sums=sums+(1-prm.gama)*cellCosts{prm.m};
% This should be the cost of context and content in single views
sums=0;
for i=1:prm.m
sums=sums+Weights(i)*cellCosts{i};
end
[co,new_t]=min(sums);
end
