%% COMMUNITIES IN THE LAST 10 YEARS

%% addpath
addpath('/Users/marijamaneva/Desktop/work/BCT/2019_03_03_BCT')

%% anni
anni = 1970:2021;
%% container
Modularity = zeros(length(anni), 1);
NumComm = zeros(length(anni), 1);
Coreness = zeros(length(anni), 1);
% Add the following container for community information
CommunityDensityInfo = cell(length(anni), 2);

%% main

for i = 1:length(anni)
    eval(['Adj = adjacency_matrix_attack_', num2str(anni(i)), ';']); %contiene semplicemente tutte le matrici e sceglie quella specificata da i
    eval(['Nodi = nodes_attack_',num2str(anni(i)),';']); %uguale per le matrici --> lo puoi togliere
    tolgo = intersect(find(sum(Adj) == 0), find(sum(Adj, 2) == 0));
    Adj(tolgo, :) = [];
    Adj(:, tolgo) = [];
    
    if ~isempty(Adj) 

        [Ci, Q] = modularity_dir(Adj); % compute community detection 
        Modularity(i) = Q; %modularity
        NumComm(i) = length(unique(Ci)); %number of communities
        
        %---------------------------------------------------------------%
        % NODES AND THEIR DENSITIES 
        % Compute community density and store both density and names
        numNodes = length(Ci);
        CommunityDensity = zeros(1, NumComm(i));
        CommunityNames = cell(1, NumComm(i));

        for j = 1:NumComm(i)
            communityNodes = find(Ci == j);
            communityAdj = Adj(communityNodes, communityNodes);
            numEdges = nnz(communityAdj) / 2; % divide by 2 to get undirected edges
            
            % Check if there are nodes in the community
            if numNodes > 1
                maxPossibleEdges = (numNodes * (numNodes - 1)) / 2;
                CommunityDensity(j) = numEdges / maxPossibleEdges;
            else
                % Handle the case where there are no or only one node in the community
                CommunityDensity(j) = 0;
            end
            
            CommunityNames{j} = strcat('Community_', num2str(j));
        end

        CommunityDensityInfo{i, 1} = CommunityDensity;
        CommunityDensityInfo{i, 2} = CommunityNames;
        

        %---------------------------------------------------------------%
        % Container for storing information about each community 
        % NAMES AND COMMUNITIES
        CommunityNodesInfo = cell(NumComm(i), 1);

        % Iterate over communities and get names of nodes
        for j = 1:NumComm(i)
            communityNodes = find(Ci == j);
            CommunityNodesInfo{j} = Nodi(communityNodes);
        end

        CommunityInfo{i} = CommunityNodesInfo;
        
        %---------------------------------------------------------------%        
        % VISUAL INSPECTION OF COMMUNITIES IN YEARS 2008,2014,2015
        % visually inspect network partition (only for 1 time stamp)
        if i==39 
            [On,Wr] = reorder_mod(Adj,Ci);
                G1 = graph(Adj);
                G1.Nodes.Name = nodes_attack_2008; %it takes the year it is referring to with i
                p1=plot(G1,'layout','force');
                title('Communities in 2008');
                p1.NodeCData = Ci;        
                axis off
                p1.MarkerSize=8;                                
        end  
        
        hold on;  
        
        if i==45
            [On,Wr] = reorder_mod(Adj,Ci);
                G2 = graph(Adj);
                G2.Nodes.Name = nodes_attack_2014; %it takes the year it is referring to with i
                figure
                p2=plot(G2,'layout','force');
                title('Communities in 2014');
                p2.NodeCData = Ci;        
                axis off
                p2.MarkerSize=8;                                
        end  
        
        hold on;
        
         if i==46
            [On,Wr] = reorder_mod(Adj,Ci);
                G3 = graph(Adj);
                G3.Nodes.Name = nodes_attack_2015; %it takes the year it is referring to with i
                figure
                p3=plot(G3,'layout','force');
                title('Communities in 2015');
                p3.NodeCData = Ci;        
                axis off
                p3.MarkerSize=8;                                
         end  
        
         hold on;
         
         if i==52
            [On,Wr] = reorder_mod(Adj,Ci);
                G4 = graph(Adj);
                G4.Nodes.Name = nodes_attack_2021; %it takes the year it is referring to with i
                figure
                p4=plot(G4,'layout','force');
                title('Communities in 2021');
                p4.NodeCData = Ci;        
                axis off
                p4.MarkerSize=10;                                
         end  
           
    end
end



%% COMMUNITIES PLOTS
figure()
subplot(2,1,1)
 
plot(Modularity, 'b')
hold on
plot(movmean(Modularity, 1), 'r')
axis tight
grid on
title('Modularity-Years moving average')
ylabel('Modularity')
 
% Adjust x-axis ticks and labels
xticks(1:length(Modularity))
xticklabels(anni) % Use the entire range of years
xtickangle(30)
 
 
% another plot
subplot(2,1,2)
 
plot(Modularity, 'b')
hold on
plot(movmean(Modularity, 1), 'r')
axis tight
grid on
title('Modularity-5 Years moving average')
ylabel('Modularity')
 
% Adjust x-axis ticks and labels for every 5 years
yearsToShow = 1970:5:2021;
xticks(1:5:length(Modularity))
xticklabels(yearsToShow)
xtickangle(30)
 
 
figure() % subplot makes more plot in one figure see the subplot function in matlab help
plot(NumComm,'b') % see plot, axis and grid functions in help
hold on
plot(movmean(NumComm,1),'r') %plot year moving average of the quantity
axis tight
grid on
title('Number of Community-Years moving average')
ylabel('Number of Community')
% Adjust x-axis ticks and labels
xticks(1:length(NumComm))
xticklabels(anni) % Use the entire range of years
xtickangle(30) % rotate the label for avoiding overlal
 
%% NODES IN THE MOST DENSE COMMUNITIES IN 2008, 2014, 2015
% Years of interest
yearsOfInterest = [2008, 2014, 2015];

% Find and display communities with highest densities for the specified years
for i = 1:length(yearsOfInterest)
    yearIndex = find(anni == yearsOfInterest(i));

    if ~isempty(CommunityDensityInfo{yearIndex, 1})
        densities = CommunityDensityInfo{yearIndex, 1};
        names = CommunityDensityInfo{yearIndex, 2};

        % Find the community with the highest density
        [~, maxDensityIndex] = max(densities);

        % Display information for the community with the highest density
        fprintf('Year: %d\n', yearsOfInterest(i));
        fprintf('Community with the highest density:\n');
        fprintf('Community Name: %s\n', names{maxDensityIndex});

        % Display node names in the community
        fprintf('Node Names:\n');
        communityNodes = find(Ci == maxDensityIndex);
        for nodeIndex = 1:length(communityNodes)
            fprintf('%s\n', Nodi{communityNodes(nodeIndex)});
        end

        fprintf('\n');
    else
        fprintf('No communities found for year %d\n\n', yearsOfInterest(i));
    end
end

%% NODES IN THE MOST DENSE COMMUNITY BETWEEN THE 3 YEARS
% Years of interest
yearsOfInterest = [2008, 2014, 2015];

% Initialize variables to store information about the most dense community
maxDensity = 0;
maxDensityCommunity = [];
MaxDensityNodes = {};  % Variable to store node names
maxDensityYear = 0;

% Find and display the most dense community among the specified years
for i = 1:length(yearsOfInterest)
    yearIndex = find(anni == yearsOfInterest(i));

    if ~isempty(CommunityDensityInfo{yearIndex, 1})
        densities = CommunityDensityInfo{yearIndex, 1};
        names = CommunityDensityInfo{yearIndex, 2};

        % Find the community with the highest density for the current year
        [currentMaxDensity, maxDensityIndex] = max(densities);

        % Update information if the current community has a higher density
        if currentMaxDensity > maxDensity
            maxDensity = currentMaxDensity;
            maxDensityCommunity = names{maxDensityIndex};

            % Get node names in the community
            communityNodes = find(Ci == maxDensityIndex);
            MaxDensityNodes = Nodi(communityNodes);  % Save node names
            
            % Update the year
            maxDensityYear = yearsOfInterest(i);
        end
    end
end

% Display information for the most dense community
fprintf('Most Dense Community:\n');
fprintf('Community Name: %s\n', maxDensityCommunity);

% Display node names in the most dense community
fprintf('Node Names:\n');
for nodeIndex = 1:length(MaxDensityNodes)
    fprintf('%s\n', MaxDensityNodes{nodeIndex});
end

% Indicate the year
fprintf('Year: %d\n', maxDensityYear);

save('DensityNodes.mat','MaxDensityNodes');


%% NODES IN THE LAST 20 YEARS
% anni
anni2 = 2002:2021; % Adjust the years for the last 20 years

% container
Modularity2 = zeros(length(anni2), 1);
NumComm2 = zeros(length(anni2), 1);

% Container for storing community information
CommunityInfo = cell(length(anni2), 1);

%% main

for i = 1:length(anni2)
    eval(['Adj = adjacency_matrix_attack_', num2str(anni2(i)), ';']); % contiene semplicemente tutte le matrici e sceglie quella specificata da i
    eval(['Nodi = nodes_attack_', num2str(anni2(i)), ';']); % uguale per le matrici --> lo puoi togliere
    tolgo = intersect(find(sum(Adj) == 0), find(sum(Adj, 2) == 0));
    Adj(tolgo, :) = [];
    Adj(:, tolgo) = [];

    if ~isempty(Adj)
        [Ci, Q] = modularity_dir(Adj); % compute community detection 
        Modularity2(i) = Q; % modularity
        NumComm2(i) = length(unique(Ci)); % number of communities

        % Container for storing information about each community
        CommunityNodesInfo = cell(NumComm2(i), 1);

        % Iterate over communities and get names of nodes
        for j = 1:NumComm(i)
            communityNodes = find(Ci == j);
            CommunityNodesInfo{j} = Nodi(communityNodes);
        end

        CommunityInfo{i} = CommunityNodesInfo;
    end
end

% Now, CommunityInfo contains the names of nodes for each community for the last 20 years.
% Access it like CommunityInfo{year}{communityIndex} for the names of nodes in a specific community.
save('CommunityNodes20Years.mat', 'CommunityInfo');