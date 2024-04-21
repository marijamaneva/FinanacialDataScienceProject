
clc; clearvars; close all;
data = load('yearly_attack_adj.mat');
attacks_per_year = data.adjacency_attacks;

figure;

% Get the number of targets
num_targets = size(attacks_per_year{1}, 2);

% Initialize a matrix to store the total attacks per target over the years
total_attacks_per_target = zeros(num_targets, length(attacks_per_year));

% Iterate through each year
for year = 1:length(attacks_per_year)
    if year == 24
        continue;
    end 
    % Get the matrix for the current year
    current_year_matrix = attacks_per_year{year};
    
    % Sum the number of attacks per target (sum along rows)
    attacks_per_target = sum(current_year_matrix, 1);
    
    % Store the results for the current year
    total_attacks_per_target(:, year) = attacks_per_target';
end
total_attacks_per_target( :,24) = [];

% Create a bar plot with different colors
b= bar(total_attacks_per_target', 'stacked');
colors = hsv(length(data.targets)); % You can use any colormap you prefer
% Assign colors to each set of bars representing a target type
for k = 1:length(data.targets)
    set(b(k),'FaceColor', colors(k,:));
end


% Add labels and legend
xlabel('Year');
ylabel('Number of Attacks');
title('Number of Attacks per Year and Target');

years = 1970:2021;
years(:, any(years == 1993, 1)) = [];
xticks(1:51);
xticklabels(years);  % Replace with the actual range of years in your data
legend(data.targets); % Add legends for each target

% Display the plot

