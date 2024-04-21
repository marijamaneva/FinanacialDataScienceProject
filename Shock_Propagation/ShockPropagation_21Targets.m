close all
clear variables
clc

%% DATA
load('TerroristProjection_20years_21Targets.mat')
load('TerroristGroups_20years.mat')

%% Shock Propagation
Adj = terrorist_projection;
iter=100;
Shockmat=zeros(size(Adj,1), iter);
initshock=ismember(terroristGroups,["Al-Shabaab","Boko Haram"]);
Shockmat(initshock,1)=1;
P=Adj./repmat(sum(Adj,2),1,size(Adj,2));
for i=1:iter-1
    Shockmat(:,i+1)=P*Shockmat(:,i);
    Shockmat(initshock,i+1)=Shockmat(initshock,1);
end
Shockmat(initshock,:)=[];

% Plot influence levels of different nodes/groups over time
figure;
plot(Shockmat');
title('Influence Levels Over Time');
xlabel('Time Steps');
ylabel('Influence');

%% Community Detection

% Perform community detection using Louvain algorithm
community = community_louvain(Adj);

%% Creating each community with names 
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

%% Community analysis

G = graph(Adj);

degreeCentrality = centrality(G, 'degree');
betweennessCentrality = centrality(G, 'betweenness');

% Identify highly central nodes and their communities
threshold = 0.9 * max(degreeCentrality); 
highlyCentralNodes = find(degreeCentrality > threshold); 

% Find communities of highly central nodes
communitiesOfCentralNodes = community(highlyCentralNodes);

%% Shock Propagation
Adj = terrorist_projection; % load terrorist interactions
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
Shockmat(initshock,:)=[];
plot(Shockmat');
title('Influence Levels Over Time');
xlabel('Time Steps');
ylabel('Influence');

% %%
% % Initialize shock matrix for each community
% ShockmatByCommunity = zeros(numCommunities, iter);
% for i = 1 : numCommunities
%     idx = find(community == i);
%     % Aggregate shock propagation values for the current community
%     ShockmatByCommunity(i, :) = mean(Shockmat(idx,:));
% end
% 
% % Plot influence levels of different nodes/groups over time
% figure;
% plot(ShockmatByCommunity');
% title('Influence Levels Over Time');
% xlabel('Time Steps');
% ylabel('Influence');
% 
