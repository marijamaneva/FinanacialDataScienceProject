clc; clearvars; close all;

%% load data
network = load('yearly_attack_adj.mat');
adjacency_attacks = network.adjacency_attacks;
terrorists = network.terrorists;
t_deg_cent=[];
t_bet_cent=[];
t_eig_cent=[];
top_deg = zeros(187, 1);
top_bet = zeros(187, 1);
top_eig = zeros(187, 1);

top_deg_l20 = zeros(187, 1);
top_bet_l20 = zeros(187, 1);
top_eig_l20 = zeros(187, 1);

%% calculate centralities
for i=0:51
    if i==23
        continue;
    end
    g = graph(adjacency_attacks{i+1}*adjacency_attacks{i+1}');
   
    % Calculate degree centrality
    degreeCentrality = centrality(g, 'degree');
    t_deg_cent= [t_deg_cent,degreeCentrality];

    % Calculate betweenness centrality
    betweennessCentrality = centrality(g, 'betweenness');
    t_bet_cent=[t_bet_cent,betweennessCentrality];
    
    % Calculate closeness centrality
    eigenCentrality = centrality(g, 'eigenvector');
    t_eig_cent=[t_eig_cent,eigenCentrality];
   
end

%% last 20 years calculation
for i=31:51
    % degree cenytrality 
    degreeCentrality = t_deg_cent(:,i);

    % Find the indices of the top 40 centralities
    [~, topIndices] = sort(degreeCentrality, 'descend');
    topIndices = topIndices(1:40);
    
    % Update countTop10 array for the top 40 indices
    top_deg_l20(topIndices) = top_deg_l20(topIndices) + 1;

    % betweenness cenytrality 
    betweennessCentrality = t_bet_cent(:,i);

     % Find the indices of the top 40 centralities
    [~, topIndices] = sort(betweennessCentrality, 'descend');
    topIndices = topIndices(1:40);
    
    % Update countTop10 array for the top 40 indices
    top_bet_l20(topIndices) = top_bet_l20(topIndices) + 1;

    % eigen cenytrality 
    eigenCentrality = t_eig_cent(:,i);

     % Find the indices of the top 40 centralities
    [~, topIndices] = sort(eigenCentrality, 'descend');
    topIndices = topIndices(1:40);
    
    % Update countTop10 array for the top 40 indices
    top_eig_l20(topIndices) = top_eig_l20(topIndices) + 1;

end


% Choose 40 most appearing terrorists
[~, sortedIndices] = sort(top_deg_l20, 'descend');
top40IndicesDeg = sortedIndices(1:40);


% Choose 40 most appearing terrorists
[~, sortedIndices] = sort(top_bet_l20, 'descend');
top40IndicesBet = sortedIndices(1:40);

% Choose 40 most appearing terrorists
[~, sortedIndices] = sort(top_eig_l20, 'descend');
top40IndicesEig = sortedIndices(1:40);

unionValues = union(union(top40IndicesDeg,top40IndicesBet,'rows','legacy'),top40IndicesEig,'rows','legacy');


%% save adj matrices
save('top_influential_terrorist.mat', 'top40IndicesDeg','top40IndicesBet','top40IndicesEig','unionValues');