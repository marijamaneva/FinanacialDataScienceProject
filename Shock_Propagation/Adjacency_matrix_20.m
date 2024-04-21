% In this script we generate 1 terrorist-target network (21 target type)
% and 1 terrorist-target_nationalities network for the last 20 years with
% their respective terrorist-terrorist projections.

close all
clear variables
clc

%% DATA FOR YEARS 2001-2021
data = readtable('processedGTD.csv');
last_20 = find(data.iyear > 2000);
data = data(last_20,:);


%% Adjacency Matrix with 21 Target Types
terroristGroups = unique(data.gname);
targetTypes = unique(data.targtype1_txt);
numGroups = numel(terroristGroups);
numTargets = numel(targetTypes);

% Create adjacency matrix for the current year    
adjacencyMatrix = zeros(numGroups, numTargets);

% Loop through the dataset for the current year
for i = 1:numGroups
    group = terroristGroups{i};
    groupAttacks = data(strcmp(data.gname, group), :);

    for j = 1:numTargets
        targetType = targetTypes{j};
        attacksOnTarget = groupAttacks(strcmp(groupAttacks.targtype1_txt, targetType), :);

        % Calculate the weighted directed edges based on the number of attacks
        numAttacks = numel(attacksOnTarget); % Weighted edge count
        
        % Update the adjacency matrix
        adjacencyMatrix(i, j) = numAttacks;
    end
end

%% Projections 21 target Types

Adj = adjacencyMatrix;

% Create target-target projection
terrorist_projection = Adj * Adj'; % Projecting targets by multiplying with its transpose

% Create terrorist-terrorist projection
target_projection = Adj' * Adj; % Projecting terrorists by multiplying with its transpose

save("TerroristProjection_20years_21Targets.mat","terrorist_projection");
save("TerroristGroups_20years.mat","terroristGroups");



%% Adjacency Matrix - Target Types combined with Nationalities
terroristGroups = unique(data.gname);
targetTypes = unique(data.combined);
numGroups = numel(terroristGroups);
numTargets = numel(targetTypes);

% Create adjacency matrix for the current year    
adjacencyMatrix = zeros(numGroups, numTargets);

% Loop through the dataset for the current year
for i = 1:numGroups
    group = terroristGroups{i};
    groupAttacks = data(strcmp(data.gname, group), :);

    for j = 1:numTargets
        targetType = targetTypes{j};
        attacksOnTarget = groupAttacks(strcmp(groupAttacks.combined, targetType), :);

        % Calculate the weighted directed edges based on the number of attacks
        numAttacks = numel(attacksOnTarget); % Weighted edge count
        
        % Update the adjacency matrix
        adjacencyMatrix(i, j) = numAttacks;
    end
end

%% Projection - Terrorist Projection With Target Types combined with Nationalities

Adj = adjacencyMatrix;

% Create target-target projection
terrorist_projection = Adj * Adj'; % Projecting targets by multiplying with its transpose

% Create terrorist-terrorist projection
target_projection = Adj' * Adj; % Projecting terrorists by multiplying with its transpose

save("terrorists_projection_20.mat","terrorist_projection");
