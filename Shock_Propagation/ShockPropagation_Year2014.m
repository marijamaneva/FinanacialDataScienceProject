% In this script we retrieve the terrorist network for a specific year
% (2014) from the matrices generated in the 'Projection' script. After
% a quick visualization using the graph command, we perform basic centrality
% and community analyses. Upon the results, we select the most influencial
% terrorist group and community to perform shock propagation. Then, we find
% which groups are more influenced when the influencial ones get shocked. 

clc
clear variables
close all
%% load similarity
load terrorists_years.mat
load terrorists_projection_yearly.mat

%% Choose the terrorist group or community to shock (e.g., first group)
YearIndex = 44; % year 2014
Adj = terrorist_proj{YearIndex}; 
terroristGroups = terroristGroups_years{YearIndex,1};

%% Terrorist-Terrorist network for Selected Year

% Set the diagonal elements to zero
Adj = Adj - diag(diag(Adj));

% Create a graph object
G = graph(Adj);

% Assign node names to the graph
G.Nodes.Name = terroristGroups;

% Find nodes with zero degree (no outgoing/ingoing interactions)
zeroDegreeNodes = find(degree(G) == 0);

% Remove nodes with zero outdegree
G = rmnode(G, zeroDegreeNodes);
terroristGroups(zeroDegreeNodes) = [];
Adj(zeroDegreeNodes,:) = [];
Adj(:,zeroDegreeNodes) = [];

% Get the degrees of all nodes
nodeDegrees = degree(G);

% Find the indices of the top 5 nodes with the highest degree
[sortedDegrees, sortedIndices] = sort(nodeDegrees, 'descend');
top5Nodes = sortedIndices(1:5);
otherNodes = sortedIndices(6:end);

% % Calculate node sizes based on degrees
% minNodeSize = 4;
% maxNodeSize = 10;
% scaledNodeSizes = (nodeDegrees - min(nodeDegrees)) / (max(nodeDegrees) - min(nodeDegrees)) * (maxNodeSize - minNodeSize) + minNodeSize;

% Update the plot with the modified graph
h = plot(G, 'Layout', 'force');
%h.NodeLabel = {};
h.NodeFontSize = 10;
h.NodeLabelColor = 'r';
labelnode(h,otherNodes,'')
highlight(h, top5Nodes, 'NodeColor', 'red');
% h.MarkerSize = scaledNodeSizes;
% nodeSizes = 5 + 10 * (nodeDegrees - min(nodeDegrees)) / (max(nodeDegrees) - min(nodeDegrees)); 
% h.NodeCData = nodeSizes;
title('Terrorist Interactions ')

%% Centrality Analysis

% Degree Centrality: Nodes with higher degree centrality are more connected within the network.
degreeCentrality = centrality(G, 'degree');

% Betweenness Centrality: Nodes with higher betweenness centrality act as crucial bridges in the network.
betweennessCentrality = centrality(G, 'betweenness');

% Closeness Centrality: Nodes with higher closeness centrality are closer to other nodes in the network.
closenessCentrality = centrality(G, 'closeness');

% Eigenvector Centrality: Nodes with higher eigenvector centrality are connected to other highly central nodes.
eigenvectorCentrality = centrality(G, 'eigenvector');

% Centrality Visualization
figure;

% Define the colors for bars
baseColor = [0.3 0.5 0.7]; % Base color for bars
highlightColor = [1 0 0];   % Color for highlighting the maximum value bar

% Degree Centrality subplot
subplot(2, 2, 1);
bar(degreeCentrality, 'FaceColor', baseColor);
hold on;
[val, idx] = max(degreeCentrality);
nameHighestDegree = terroristGroups{idx};
bar(idx, val, 'FaceColor', highlightColor);
ylabel('Degree Centrality');
xlabel('Terrorist Groups');
xticks(idx);
xticklabels(terroristGroups(idx));
fprintf('Terrorist Group with Highest Degree Centrality: %s (Centrality: %.4f)\n', nameHighestDegree, val);

% Betweenness Centrality subplot
subplot(2, 2, 2);
bar(betweennessCentrality, 'FaceColor', baseColor);
hold on;
[val, idx] = max(betweennessCentrality);
nameHighestBetweenness = terroristGroups{idx};
bar(idx, val, 'FaceColor', highlightColor);
ylabel('Betweenness Centrality');
xlabel('Terrorist Groups');
xticks(idx);
xticklabels(terroristGroups(idx));
fprintf('Terrorist Group with Highest Betweenness Centrality: %s (Centrality: %.4f)\n', nameHighestBetweenness, val);


% Closeness Centrality subplot
subplot(2, 2, 3);
bar(closenessCentrality, 'FaceColor', baseColor);
hold on;
[val, idx] = max(closenessCentrality);
nameHighestCloseness = terroristGroups{idx};
bar(idx, val, 'FaceColor', highlightColor);
ylabel('Closeness Centrality');
xlabel('Terrorist Groups');
xticks(idx);
xticklabels(terroristGroups(idx));
fprintf('Terrorist Group with Highest Closeness Centrality: %s (Centrality: %.4f)\n', nameHighestCloseness, val);

% Eigenvector Centrality subplot
subplot(2, 2, 4);
bar(eigenvectorCentrality, 'FaceColor', baseColor);
hold on;
[val, idx] = max(eigenvectorCentrality);
nameHighestEigenvector = terroristGroups{idx};
bar(idx, val, 'FaceColor', highlightColor);
ylabel('Eigenvector Centrality');
xlabel('Terrorist Groups');
xticks(idx);
xticklabels(terroristGroups(idx));
fprintf('Terrorist Group with Highest Eigenvector Centrality: %s (Centrality: %.4f)\n', nameHighestEigenvector, val);



%% Community Analysis

% Perform community detection using Louvain algorithm
community = community_louvain(Adj);

% Get the number of nodes and communities
numNodes = numnodes(G);
numCommunities = max(community);

% Calculate community sizes (number of nodes in each community)
communitySizes = zeros(1, numCommunities);
for i = 1:numCommunities
    communitySizes(i) = sum(community == i);
end

% Calculate community densities (average edge density within each community)

communityDensities = zeros(1, numCommunities);
for i = 1:numCommunities
    nodesInCommunity = find(community == i);
    subgraph = Adj(nodesInCommunity, nodesInCommunity);
    numEdges = nnz(subgraph) / 2; % Divide by 2 as the graph is undirected
    numPossibleEdges = numel(subgraph) / 2; % Total possible edges in an undirected graph
    communityDensities(i) = numEdges / numPossibleEdges;
end


centralNodes = zeros(1, numCommunities);
for i = 1:numCommunities
    nodesInCommunity = find(community == i);
    subgraph = Adj(nodesInCommunity, nodesInCommunity);
    
    % Calculate degree centrality for nodes in the community subgraph
    subgraphG = graph(subgraph);
    degreeCentrality2 = centrality(subgraphG, 'degree');
    
    [~, idx] = max(degreeCentrality2);
    centralNodes(i) = nodesInCommunity(idx);
end


% Index of the largest community's central node
idxLargestCommunityCentral = find(communitySizes == max(communitySizes), 1);
nameLargestCommunityCentral = terroristGroups{centralNodes(idxLargestCommunityCentral)};

% Index of the community with highest density's central node
idxHighestDensityCentral = find(communityDensities == max(communityDensities), 1);
nameHighestDensityCentral = terroristGroups{centralNodes(idxHighestDensityCentral)};

% Display names associated with central nodes
fprintf('Terrorist Group associated with central node of the largest community: %s\n', nameLargestCommunityCentral);
fprintf('Terrorist Group associated with central node of the community with highest density: %s\n', nameHighestDensityCentral);

% Create a cell structure to store groups in each community
groups_in_community = cell(max(community), 1);
numCommunities = max(community);
% Map groups to their respective communities
for idx = 1:length(terroristGroups)
    comm = community(idx);
    if isempty(groups_in_community{comm})
        groups_in_community{comm} = terroristGroups(idx);
    else
        groups_in_community{comm} = [groups_in_community{comm}; terroristGroups(idx)];
    end
end

% Identify highly central nodes and their communities
threshold = 0.75 * max(degreeCentrality); 
highlyCentralNodes = find(degreeCentrality > threshold); 
disp('Highly central nodes:');
disp(terroristGroups(highlyCentralNodes));

% Find communities of highly central nodes
communitiesOfCentralNodes = community(highlyCentralNodes);
disp('Communities of highly central nodes:');
disp(communitiesOfCentralNodes);

% % Find communities of the 5 most highly central nodes
% disp('Top 5 highly central nodes:');
% disp(terroristGroups(top5Nodes));
% communitiesOfCentralNodes = community(top5Nodes);
% disp('Communities of top 5 highly central nodes:');
% disp(communitiesOfCentralNodes);

%% Simulation of shock propagation - highest degree centrality
Adj = terrorist_proj{YearIndex}; 
terroristGroups = terroristGroups_years{YearIndex,1};
iter=100;
Shockmat=zeros(size(Adj,1), iter);
initshock=ismember(terroristGroups_years{YearIndex},nameHighestDegree);
Shockmat(initshock,1)= 1;
P=Adj./repmat(sum(Adj,2),1,size(Adj,2));
for i=1:iter-1
    Shockmat(:,i+1)=P*Shockmat(:,i);
    Shockmat(initshock,i+1)=Shockmat(initshock,1);
end
Shockmat(initshock,:)=[];
terroristGroups(initshock) = [];

% Plot influence levels of different nodes/groups over time
figure;
plot(Shockmat');
title('Influence Levels Over Time');
xlabel('Time Steps');
ylabel('Influence');

%% Identify Most Highly Impacted Groups
finalShock = Shockmat(:, end); % Extract shock state at the final iteration
[sortedShock, groupIndices] = sort(finalShock, 'descend'); % Sort shock values in descending order

% Display the top three most highly impacted groups
over90 = length(find(sortedShock > 0.90))
indices_of_highest_ones = groupIndices(1:over90);
terrorists = terroristGroups;
% terrorists(initshock) = [];
highlyImpactedGroups = terrorists(indices_of_highest_ones)% Extract top three groups
topShockValues = sortedShock(1:over90) % Corresponding shock values
meanShockValues = mean(Shockmat, 2); % Calculate mean shock value for each group
figure;
bar(meanShockValues);
title('Mean Shock Values across Terrorist Groups');
xlabel('Terrorist Groups');
ylabel('Mean Shock Values');

% Filter terrorist groups with mean shock values greater than 0.5
selectedIndices = sort(indices_of_highest_ones);
selectedLabels = terrorists(selectedIndices);

% Set x-axis tick labels only for selected terrorists
xticks(selectedIndices);
xticklabels(selectedLabels);
xtickangle(45)
hold on
bar(selectedIndices,meanShockValues(selectedIndices), 'FaceColor', [1 0 0]);



%% Simulation of shock propagation - largest community
terroristGroups = terroristGroups_years{YearIndex,1};
iter=100;
Shockmat=zeros(size(Adj,1), iter);
initiallyInfectedCommunities=mode(communitiesOfCentralNodes);
communityNodes = find(community == initiallyInfectedCommunities);
initshock = communityNodes;
Shockmat(initshock,1)=1;
P=Adj./repmat(sum(Adj,2),1,size(Adj,2));
for i=1:iter-1
    Shockmat(:,i+1)=P*Shockmat(:,i);
    Shockmat(initshock,i+1)=Shockmat(initshock,1);
end

% Plot influence levels of different nodes/groups over time
figure;
plot(Shockmat');
title('Shock Propagation - community with highly central nodes');
xlabel('Time Steps');
ylabel('Influence');

% finalShock = Shockmat(:, end); % Extract shock state at the final iteration
% finalShock(initshock) = [];
% terroristGroups(initshock) = [];
% [sortedShock, groupIndices] = sort(finalShock, 'descend'); % Sort shock values in descending order
% 
% % Display the top three most highly impacted groups
% over90 = length(find(sortedShock > 0.90));
% indices_of_highest_ones = groupIndices(1:over90);
% terrorists = terroristGroups;
% highlyImpactedGroups = terrorists(indices_of_highest_ones)


%% Shock Propagation by Community
% Initialize shock matrix for each community
ShockmatByCommunity = zeros(numCommunities, iter);
for i = 1 : numCommunities
    idx = find(community == i);
    % Aggregate shock propagation values for the current community
    ShockmatByCommunity(i, :) = mean(Shockmat(idx,:));
end

% Plot influence levels of different nodes/groups over time
figure;
plot(ShockmatByCommunity');
title('Shock Propagation - Communities');
xlabel('Time Steps');
ylabel('Influence');


%% Identify Most Highly Impacted Communities
% finalShock_com = ShockmatByCommunity(:, end); % Extract shock state at the final iteration
% [sortedShock, groupIndices] = sort(finalShock_com, 'descend'); % Sort shock values in descending order
% 
% % Display the top three most highly impacted groups
% over90 = length(find(sortedShock > 0.90))
% indices_of_highest_ones = groupIndices(1:over90);
% for i = 1:length(indices_of_highest_ones)
%     highlyImpactedCommunities = groups_in_community{indices_of_highest_ones(i)}% Extract top three groups
%     disp(highlyImpactedCommunities)
% end
