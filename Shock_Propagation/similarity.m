clc
clear variables
close all
%% load similarity
load Adjacency_matrices.mat
load terrorists_years.mat
load targets_years.mat
load terrorists_proj.mat
%% Choose the terrorist group or community to shock (e.g., first group)
YearIndex = 44;
adjacency_matrix = terrorist_proj{YearIndex}; % Retrieve the adjacency matrix for the selected year
terroristGroups = terroristGroups_years{YearIndex,1};
%%
% Calculate the similarity matrix based on common neighbors
num_nodes = size(adjacency_matrix, 1);
similarity_matrix = zeros(num_nodes);

for i = 1:num_nodes
    for j = 1:num_nodes
        if i ~= j
            common_neighbors = sum(adjacency_matrix(i, :) & adjacency_matrix(j, :));
            similarity_matrix(i, j) = common_neighbors;
        end
    end
end

%%
Adj = similarity_matrix;
iter=100;
Shockmat=zeros(size(Adj,1), iter);
initshock=contains(terroristGroups,["Muslim extremists"]);
Shockmat(initshock,1)=1;
P=Adj./repmat(sum(Adj,2),1,size(Adj,2));
for i=1:iter-1
    Shockmat(:,i+1)=P*Shockmat(:,i);
    Shockmat(initshock,i+1)=Shockmat(initshock,1);
end
% Shockmat(initshock,:)=[];

% Plot influence levels of different nodes/groups over time
figure;
plot(Shockmat');
title('Influence Levels Over Time');
xlabel('Time Steps');
ylabel('Influence');
legend(terroristGroups,"Location","bestoutside"); 

