clc; clearvars; close all;

%% Load dataset
data = readtable('processedGTD.csv');
adjacency_kills = cell(length(unique(data.iyear)),1);
target_nodes = cell(length(unique(data.iyear)),1);
terrorist_nodes = cell(length(unique(data.iyear)),1);
terrorists = unique(data.gname);
targets = unique(data.targtype1_txt);

%% Populate the adjacency matrix with the number of attacks, each year

for j=0:length(unique(data.iyear))
    if j==23
        continue;
    end
    condition = (data.iyear == (1970 + j));
    subdata = data(condition, :);
    
    % Create an adjacency matrix for terrorists and targets
    adjacencyMatrix = zeros(length(terrorists), length(targets));
    for i = 1:size(subdata, 1)
        terroristIndex = find(strcmp(subdata.gname{i}, terrorists));
        targetIndex = find(strcmp(subdata.targtype1_txt{i}, targets));
        adjacencyMatrix(terroristIndex, targetIndex) = adjacencyMatrix(terroristIndex, targetIndex) + subdata.nkill(i) ;
    end
    target_nodes{j+1} = unique(subdata.targtype1_txt);
    terrorist_nodes{j+1} = unique(subdata.gname);
    adjacencyMatrix(isnan(adjacencyMatrix)) = 0;
    adjacency_kills{j+1} = adjacencyMatrix;

end 


%% save adj matrices
save('yearly_kill_adj.mat', 'adjacency_kills','terrorist_nodes','target_nodes','terrorists');








