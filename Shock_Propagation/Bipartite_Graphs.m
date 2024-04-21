%% BIPARTITE NETWORK - BIPARTITE GRAPHS

% In this script we generate bipartite networks for each year from 1970 to
% 2021 and save the matrices as well as the terrorist groups and target
% types combined with their nationalities for each year.

% In bipartite graphs, nodes are separated into two disjoint
% sets such that links only connect nodes in different partitions. In our
% case, the primary set of nodes represents the terrorist groups, while
% the secondary set refers to the targets which constitute the objectives
% of the attacks. In particular, the targets have been combined with their
% nationalities to have better understanding.
% We generate the terrorist-target bipartite network where links represent
% the total number of yearly attacks perpetrated by a terrorist 𝑣 on a
% certain target type 𝑢.

close all
clear variables
clc

%% DATA
data = readtable('processedGTD.csv');

%% Extract unique years
years = unique(data.iyear);

%% Create an empty cell array to store adjacency matrices for each year
adjacencyMatrices = cell(numel(years), 2);
terroristGroups_years = cell(numel(years), 2);
targets_years = cell(numel(years), 2);

%% Loop through each year and construct adjacency matrices
for y = 1:numel(years)
    yearData = data(data.iyear == years(y), :);
    
    terroristGroups = unique(yearData.gname);
    targetTypes = unique(yearData.combined);
    numGroups = numel(terroristGroups);
    numTargets = numel(targetTypes);

    % Create adjacency matrix for the current year    
    adjacencyMatrix = zeros(numGroups, numTargets);
    
    % Loop through the dataset for the current year
    for i = 1:numGroups
        group = terroristGroups{i};
        groupAttacks = yearData(strcmp(yearData.gname, group), :);

        for j = 1:numTargets
            targetType = targetTypes{j};
            attacksOnTarget = groupAttacks(strcmp(groupAttacks.combined, targetType), :);

            % Calculate the weighted directed edges based on the number of attacks
            numAttacks = numel(attacksOnTarget); % Weighted edge count
            
            % Update the adjacency matrix
            adjacencyMatrix(i, j) = numAttacks;
        end
    end
    
    % Store the adjacency matrix and the current year
    adjacencyMatrices{y, 1} = adjacencyMatrix;
    adjacencyMatrices{y, 2} = years(y);
    terroristGroups_years{y, 1} = terroristGroups;
    targets_years{y, 1} = targetTypes;
    terroristGroups_years{y, 2} = years(y);
    targets_years{y, 2} = years(y);
end


% adjacencyMatrices contains temporal adjacency matrices for each year
save("Adjacency_matrices_yearly","adjacencyMatrices");
save("terrorists_years","terroristGroups_years");
save("targets_years","targets_years");

%% Create a bipartite graph per year 

% Determine how many figures you want to group your data into
total_figures = 6;  

% Calculate how many adjacency matrices you have
num_matrices = length(adjacencyMatrices);

plots_per_figure = ceil(num_matrices / total_figures);

figure_count = 0;

for fig = 1:total_figures
    figure('Name', ['Grouped Figure ', num2str(fig)]); 
    
    for plot_num = 1:plots_per_figure
        figure_count = figure_count + 1;
        
        if figure_count > num_matrices
            break;  % Break the loop if all adjacency matrices have been displayed
        end
        
        % Extract data for the current iteration
        y = figure_count;
        adjacencyMatrix = adjacencyMatrices{y, 1};
        terroristGroups = terroristGroups_years{y, 1};
        targetTypes = targets_years{y, 1};
        
        % Find non-zero elements in the adjacency matrix
        [row, col] = find(adjacencyMatrix);
        
        % Create a directed graph (digraph) object
        G = digraph();
        
        % Add nodes for terrorist groups and target types
        G = addnode(G, numel(terroristGroups) + numel(targetTypes));
        
        % Add edges to the graph based on the adjacency matrix
        for i = 1:numel(row)
            G = addedge(G, row(i), numel(terroristGroups) + col(i));
        end
        
        % Assign node names to the graph
        G.Nodes.Name = [terroristGroups; targetTypes];
        
        % Define node colors for visualization
        node_colors = [repmat([0.7 0.2 0.2], numel(terroristGroups), 1); repmat([0.2 0.7 0.2], numel(targetTypes), 1)];
        
        subplot(sqrt(plots_per_figure), sqrt(plots_per_figure), plot_num);
        h = plot(G, 'Layout', 'force', 'NodeColor', node_colors, 'MarkerSize', [2*ones(numel(terroristGroups), 1); 2*ones(numel(targetTypes), 1)]);
        h.NodeLabel = {};
        title(['Bipartite Graph for Year ' num2str(adjacencyMatrices{y, 2})]);
    end
end

