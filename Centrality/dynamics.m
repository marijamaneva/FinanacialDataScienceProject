clc; clearvars; close all;

data = load('yearly_attack_adj.mat');

years = 1970:1970+length(data.target_nodes)-1; 
terroristsDynamics = data.terrorist_nodes; 
targetsDynamics = data.target_nodes;

terroristsDynamics(24, :) = [];
targetsDynamics(24, :) = [];
years(:, any(years == 1993, 1)) = [];
% Initialize an array to store the number of terrorists for each year
numTerroristsArray = zeros(length(years), 1);
numTargetsArray = zeros(length(years), 1);
% Loop through each cell
for year = 1:length(years)
    % Get the cell for the current year
    currentYearCellterr = terroristsDynamics{year};
    currentYearCelltar = targetsDynamics{year};
    % Count the number of terrorists for the current year
    numTerrorists = numel(currentYearCellterr);
    numTargets = numel(currentYearCelltar);
    % Store the count in the array
    numTerroristsArray(year) = numTerrorists;
    numTargetsArray(year) = numTargets;
end


% Applying a 5-year moving window
windowSize = 5;
terroristsSmoothed = movmean(numTerroristsArray, windowSize);
targetsSmoothed = movmean(numTargetsArray, windowSize);
highlight_indices = 1:length(years);
% Create a plot
figure;
yyaxis left;
plot(years, targetsSmoothed, 'Color', [0.7, 0.7, 0.7], 'LineWidth', 2);
ylabel('Targets number');
hold on;
% Overlay target on highlighted points
scatter(years(highlight_indices), targetsSmoothed(highlight_indices), 80, [0.7, 0.7, 0.7],'Marker', 'v','LineWidth', 3);

yyaxis right;
plot(years, terroristsSmoothed, 'k', 'LineWidth', 2);
ylabel('Terrorists number');

xlim([min(years), max(years)]);
xlabel('Year');
title('Dynamics of Terrorists and Targets Over Time');
hold on;
% Overlay terrorrist on highlighted points
scatter(years(highlight_indices), terroristsSmoothed(highlight_indices), 80, 'black','LineWidth', 3);

grid on;

