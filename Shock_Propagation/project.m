close all
clear variables
clc

%% DATA
gtd_dataset = readtable('gtd.csv');
data = readtable('processedGTD.csv');

%% Extract unique terrorist groups, target types, and years
terroristGroups = unique(data.gname);
targetTypes = unique(data.targtype1);
years = unique(data.iyear);

%% Load adjacencyMatrices 
% contains temporal adjacency matrices for each year
adjacencyMatrices = load("Target_Terrorist_AdjMatrix.mat");
adjacencyMatrices = struct2cell(adjacencyMatrices);
adjacencyMatrices = adjacencyMatrices{1};


%% Network Analysis
% Initialize arrays to store metrics over time
numNodes = zeros(numel(years), 1);
numEdges = zeros(numel(years), 1);
densities = zeros(numel(years), 1);

% Loop through each year's adjacency matrix
for y = 1:numel(years)
    adjacencyMatrix = adjacencyMatrices{y, 1};

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
    G.Nodes.Name = [terroristGroups; cellstr(num2str(targetTypes))];


    % Find isolated nodes (nodes without any edges)
    isolated_nodes = find(indegree(G) == 0 & outdegree(G) == 0);
    
    % Create a new graph without isolated nodes
    filtered_G = rmnode(G, isolated_nodes);

    % Calculate basic metrics
    numNodes(y) =  numnodes(filtered_G);
    numEdges(y) =  numedges(filtered_G);
    possibleEdges = numNodes(y) * (numNodes(y) - 1);
    densities(y) = numEdges(y) / possibleEdges;

end

% Plotting evolution of metrics over time
figure;
subplot(2, 1, 1);
plot(years, numNodes, 'b-o');
xlabel('Year');
ylabel('Number of Nodes');
title('Evolution of Number of Nodes');

subplot(2, 1, 2);
plot(years, numEdges, 'r-o');
xlabel('Year');
ylabel('Number of Edges');
title('Evolution of Number of Edges');

figure;
plot(years, densities, 'r-o');
xlabel('Year');
ylabel('Density');
title('Evolution of Network Density');

