function [CR, PI, T, seed] = MainIB(cellX, m, Tsize, theta, beta, Restarts, UniformPrior, CentroidFigs)
% Main Information Bottleneck (IB) algorithm implementation
% Inputs:
%   cellX - cell array containing input data
%   m - number of clusters/partitions
%   Tsize - size of the bottleneck representation
%   theta - convergence threshold parameter
%   beta - trade-off parameter between compression and relevance
%   Restarts - number of random restarts for optimization
%   UniformPrior - flag for using uniform prior distribution
%   CentroidFigs - flag for generating centroid visualization figures
%
% Outputs:
%   CR - Cluster Results matrix (Restarts x data points)
%   PI - Final optimal objective function values for each restart
%   T - Optimal bottleneck representation found
%   seed - Random seed used for the run

format compact

% Initialize algorithm parameters and process input data
[prm] = Initparameters(m, Tsize, theta, beta, Restarts, UniformPrior, CentroidFigs);
[cellInp] = ProcessInput(cellX, prm);

% Initialize result storage
CR = zeros(Restarts, size(cellX{1}, 1)); % Clustering result for each restart
PI = zeros(1, Restarts);                 % Final optimal objective value for each restart
bestL = -inf;                            % Track best objective value across all restarts

% Main optimization loop with multiple random restarts
for RS = 1:prm.Restarts
%     fprintf('Restart number %d....\n', RS);

    % Generate random initial partition and optimize
    [tmpT] = RandomPartition(cellX, cellInp, prm);
    [tmpT] = OptimizeT(tmpT, cellX, cellInp, prm);

    % Store results for current restart
    CR(RS, :) = tmpT{1}.Pt_x;      % Cluster assignments
    PI(RS) = tmpT{1}.L(end);       % Final objective value

    % Update global best solution if current restart is better
    if tmpT{1}.L(end) > bestL
        T{1} = tmpT{1};
        bestL = T{1}.L(end);       % Use final value as best (converged)
    end
end

% Optional: Generate centroid visualization figures
% (Code for centroid figures would go here if CentroidFigs is enabled)

% Calculate additional information-theoretic metrics (commented out)
% for i = 1:prm.m
%     cellT{i}.Ity_div_Ixy = cellT{i}.Ity(end) / cellInp{i}.Ixy;
%     cellT{i}.Ht_div_Hx = cellT{i}.Ht(end) / cellInp{i}.Hx;
% end

% Return the random seed used for reproducibility
seed = prm.RunSeed;

end