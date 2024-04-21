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

    % Find the indices of the top 10 centralities
    [~, topIndices] = sort(degreeCentrality, 'descend');
    topIndices = topIndices(1:10);
    
    % Update countTop10 array for the top 10 indices
    top_deg(topIndices) = top_deg(topIndices) + 1;


    % Calculate betweenness centrality
    betweennessCentrality = centrality(g, 'betweenness');
    t_bet_cent=[t_bet_cent,betweennessCentrality];
    
    % Find the indices of the top 10 centralities
    [~, topIndices] = sort(betweennessCentrality, 'descend');
    topIndices = topIndices(1:10);
    
    % Update countTop10 array for the top 10 indices
    top_bet(topIndices) = top_bet(topIndices) + 1;


    % Calculate closeness centrality
    eigenCentrality = centrality(g, 'eigenvector');
    t_eig_cent=[t_eig_cent,eigenCentrality];
   
    % Find the indices of the top 10 centralities
    [~, topIndices] = sort(eigenCentrality, 'descend');
    topIndices = topIndices(1:10);
    
    % Update countTop10 array for the top 10 indices
    top_eig(topIndices) = top_eig(topIndices) + 1;
end
%% plot top 10 degree centrality 
% Choose 10 most appearing terrorists
[~, sortedIndices] = sort(top_deg, 'descend');
top10Indices = sortedIndices(1:10);

% Extract the data for the selected rows
selectedData = t_deg_cent(top10Indices, :);

% Create a figure
figure;
xLimit = [0, 53];
yLimit = [0, 100];
% Plot each of the top 10 indices in separate subplots
for i = 1:10
    subplot(5, 2, i);  % Assuming you want a 5x2 grid of subplots
    plot(selectedData(i,:)', 'LineWidth', 2);
    % Set common axis limits
    xlim(xLimit);
    ylim(yLimit);
    title([terrorists(top10Indices(i))]);
    xlabel('year');
    ylabel('degree centrality');
    grid on;
end

% Adjust the layout
sgtitle('Top 10 terrorists hubs based on degree centrality ');
% Maximize the figure window
set(gcf, 'Position', get(0, 'Screensize'));
saveas(gcf, 'top10degcent.jpg');
%% plot top 10 betweenness centrality  
% Choose 10 most appearing terrorists
[~, sortedIndices] = sort(top_bet, 'descend');
top10Indices = sortedIndices(1:10);

% Extract the data for the selected rows
selectedData = t_bet_cent(top10Indices, :);

% Create a figure
figure;
xLimit = [0, 53];
yLimit = [0, 100];
% Plot each of the top 10 indices in separate subplots
for i = 1:10
    subplot(5, 2, i);  % Assuming you want a 5x2 grid of subplots
    plot(selectedData(i,:)', 'LineWidth', 2);
    % Set common axis limits
    xlim(xLimit);
    ylim(yLimit);
    title([terrorists(top10Indices(i))]);
    xlabel('year');
    ylabel('betweenness centrality');
    grid on;
end

% Adjust the layout
sgtitle('Top 10 terrorists hubs based on betweenness centrality ');
% Maximize the figure window
set(gcf, 'Position', get(0, 'Screensize'));
saveas(gcf, 'top10betcent.jpg');
%% plot top 10 eigen centrality  
% Choose 10 most appearing terrorists
[~, sortedIndices] = sort(top_eig, 'descend');
top10Indices = sortedIndices(1:10);

% Extract the data for the selected rows
selectedData = t_eig_cent(top10Indices, :);

% Create a figure
figure;
% Plot each of the top 10 indices in separate subplots
for i = 1:10
    subplot(5, 2, i);  % Assuming you want a 5x2 grid of subplots
    plot(selectedData(i,:)', 'LineWidth', 2);
  
    title([terrorists(top10Indices(i))]);
    xlim(xLimit);
    xlabel('year');
    ylabel('eigen centrality');
    grid on;
end

% Adjust the layout
sgtitle('Top 10 terrorists hubs based on eigen centrality ');
% Maximize the figure window
set(gcf, 'Position', get(0, 'Screensize'));
saveas(gcf, 'top10eigcent.jpg');
%% last 20 years calculation
for i=31:51
    % degree cenytrality 
    degreeCentrality = t_deg_cent(:,i);

    % Find the indices of the top 10 centralities
    [~, topIndices] = sort(degreeCentrality, 'descend');
    topIndices = topIndices(1:10);
    
    % Update countTop10 array for the top 10 indices
    top_deg_l20(topIndices) = top_deg_l20(topIndices) + 1;

    % betweenness cenytrality 
    betweennessCentrality = t_bet_cent(:,i);

     % Find the indices of the top 10 centralities
    [~, topIndices] = sort(betweennessCentrality, 'descend');
    topIndices = topIndices(1:10);
    
    % Update countTop10 array for the top 10 indices
    top_bet_l20(topIndices) = top_bet_l20(topIndices) + 1;

    % eigen cenytrality 
    eigenCentrality = t_eig_cent(:,i);

     % Find the indices of the top 10 centralities
    [~, topIndices] = sort(eigenCentrality, 'descend');
    topIndices = topIndices(1:10);
    
    % Update countTop10 array for the top 10 indices
    top_eig_l20(topIndices) = top_eig_l20(topIndices) + 1;

end
%% plot last 20 years degree centrality 

% Choose 10 most appearing terrorists
[~, sortedIndices] = sort(top_deg_l20, 'descend');
top10Indices = sortedIndices(1:10);

% Extract the data for the selected rows
selectedData = t_deg_cent(top10Indices, :);

% Create a figure
figure;
xLimit = [0, 53];
yLimit = [0, 100];
% Plot each of the top 10 indices in separate subplots
for i = 1:10
    subplot(5, 2, i);  % Assuming you want a 5x2 grid of subplots
    plot(selectedData(i,:)', 'LineWidth', 2);
    % Set common axis limits
    xlim(xLimit);
    ylim(yLimit);
    title([terrorists(top10Indices(i))]);
    xlabel('year');
    ylabel('degree centrality');
    grid on;
end

% Adjust the layout
sgtitle('Top 10 terrorists hubs in last 20 years based on degree centrality ');
% Maximize the figure window
set(gcf, 'Position', get(0, 'Screensize'));
saveas(gcf, 'top10l20degcent.jpg');
%% plot last 20 years betweenness centrality 

% Choose 10 most appearing terrorists
[~, sortedIndices] = sort(top_bet_l20, 'descend');
top10Indices = sortedIndices(1:10);

% Extract the data for the selected rows
selectedData = t_bet_cent(top10Indices, :);

% Create a figure
figure;
xLimit = [0, 53];
yLimit = [0, 100];
% Plot each of the top 10 indices in separate subplots
for i = 1:10
    subplot(5, 2, i);  % Assuming you want a 5x2 grid of subplots
    plot(selectedData(i,:)', 'LineWidth', 2);
    % Set common axis limits
    xlim(xLimit);
    ylim(yLimit);
    title([terrorists(top10Indices(i))]);
    xlabel('year');
    ylabel('betweenness centrality');
    grid on;
end

% Adjust the layout
sgtitle('Top 10 terrorists hubs in last 20 years based on betweenness centrality ');
% Maximize the figure window
set(gcf, 'Position', get(0, 'Screensize'));
saveas(gcf, 'top10l20betcent.jpg');
%% plot last 20 years eigen centrality 

% Choose 10 most appearing terrorists
[~, sortedIndices] = sort(top_eig_l20, 'descend');
top10Indices = sortedIndices(1:10);

% Extract the data for the selected rows
selectedData = t_eig_cent(top10Indices, :);

% Create a figure
figure;

% Plot each of the top 10 indices in separate subplots
for i = 1:10
    subplot(5, 2, i);  % Assuming you want a 5x2 grid of subplots
    plot(selectedData(i,:)', 'LineWidth', 2);

    title([terrorists(top10Indices(i))]);
    xlim(xLimit);
    xlabel('year');
    ylabel('eigen centrality');
    grid on;
end

% Adjust the layout
sgtitle('Top 10 terrorists hubs in last 20 years based on eigen centrality ');
% Maximize the figure window
set(gcf, 'Position', get(0, 'Screensize'));
saveas(gcf, 'top10l20eigcent.jpg');



