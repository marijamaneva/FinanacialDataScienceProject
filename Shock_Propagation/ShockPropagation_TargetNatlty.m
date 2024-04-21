close all
clear variables
clc

%% DATA
load("terrorists_projection_20.mat");
load("TerroristGroups_20years.mat");

%% Shock Propagation - Al-Shabaab and Boko Haram
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

finalShock = Shockmat(:, end); % Extract shock state at the final iteration
[sortedShock, groupIndices] = sort(finalShock, 'descend'); % Sort shock values in descending order

% Display the top three most highly impacted groups
over80 = length(find(sortedShock > 0.80))
indices_of_highest_ones = groupIndices(1:over80);
terrorists = terroristGroups;
terrorists(initshock) = [];
highlyImpactedGroups = terrorists(indices_of_highest_ones)% Extract top three groups
topShockValues = sortedShock(1:over80) % Corresponding shock values


%% Shock Propagation - Muslim Extremists and TTP
Adj = terrorist_projection;
iter=100;
Shockmat=zeros(size(Adj,1), iter);
initshock=ismember(terroristGroups,["Muslim extremists","Tehrik-i-Taliban Pakistan (TTP)"]);
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


% Identify Most Highly Impacted Groups
finalShock = Shockmat(:, end); % Extract shock state at the final iteration
[sortedShock, groupIndices] = sort(finalShock, 'descend'); % Sort shock values in descending order

% Display the top three most highly impacted groups
over80 = length(find(sortedShock > 0.80))
indices_of_highest_ones = groupIndices(1:over80);
terrorists = terroristGroups;
terrorists(initshock) = [];
highlyImpactedGroups = terrorists(indices_of_highest_ones)% Extract top three groups
topShockValues = sortedShock(1:over80) % Corresponding shock values

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