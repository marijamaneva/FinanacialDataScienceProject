% In this script we generate the terrorist-terrorist and target-target
% projections using the terrorist-target networks we obtained in the
% 'Bipartite_Graphs' script

clc
clear variables
close all
%% load similarity
load Adjacency_matrices_yearly.mat
load terrorists_years.mat
load targets_years.mat

%% Projections
years = cell2mat(adjacencyMatrices(:,2));
% Initialize variables to store projections
target_projection = cell(numel(years), 2);
terrorist_proj = cell(numel(years), 2);

for year = 1:numel(years)
    % Retrieve adjacency matrix for the current year
    Adj = adjacencyMatrices{year,1};

    % Create target-target projection
    target_projection{year,1} = Adj' * Adj; % Projecting targets by multiplying with its transpose
    target_projection{year,2} = years(year);

    % Create terrorist-terrorist projection
    terrorist_proj{year,1} = Adj * Adj'; % Projecting terrorists by multiplying with its transpose
    terrorist_proj{year,2} = years(year);
end


save("terrorists_projection_yearly","terrorist_proj");
%save("target_projection","targets_years");
