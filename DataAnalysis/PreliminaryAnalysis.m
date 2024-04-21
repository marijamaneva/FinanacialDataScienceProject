close all
clear variables
clc

%% DATA
gtd_dataset = readtable('gtd.csv');
data = readtable('processedGTD.csv');


%% Attack Location 1970-2021

% Extract latitude and longitude information
latitudes = data.latitude;
longitudes = data.longitude;
countries = data.country_txt;

% Remove missing or NaN values in latitude, longitude, and year
valid_indices = ~isnan(latitudes) & ~isnan(longitudes);
latitudes = latitudes(valid_indices);
longitudes = longitudes(valid_indices);
countries = countries(valid_indices);

tbl = tabulate(countries);
tbl = sortrows(tbl,1);
attack_counts = cell2mat(tbl(:, 2));
sorted_data = sortrows(tbl, -2);
% use only one pair of latitude and longitudes to represent a country
unique_countries = unique(countries); % Get unique country names

% Initialize arrays to store aggregated latitude and longitude data
avg_latitudes = zeros(size(unique_countries));
avg_longitudes = zeros(size(unique_countries));

% Loop through unique countries and calculate average latitude and longitude
for i = 1:numel(unique_countries)
    country = unique_countries{i};
    indices = strcmp(countries, country); % Find indices of data corresponding to the current country
    avg_latitudes(i) = mean(latitudes(indices)); % Calculate mean latitude
    avg_longitudes(i) = mean(longitudes(indices)); % Calculate mean longitude
end

% Create a new table or structure to store the unique country names along with their representative latitude and longitude
country_coordinates = table(unique_countries, avg_latitudes, avg_longitudes, 'VariableNames', {'Country', 'Latitude', 'Longitude'});

% Merge attack counts with country coordinates based on the country names
country_coordinates.attack_counts = attack_counts;
figure('Name',"Attack Location 1970-2021")
% Geographic plot
geobubble(country_coordinates.Latitude, country_coordinates.Longitude, country_coordinates.attack_counts);

% Customize the map properties

title('Attack Location 1970-2021');



%% Attack Location 2000-2021
data_recent = data(data.iyear >= 2000, :);
% Extract latitude and longitude information
latitudes = data_recent.latitude;
longitudes = data_recent.longitude;
countries = data_recent.country_txt;

% Remove missing or NaN values in latitude, longitude, and year
valid_indices = ~isnan(latitudes) & ~isnan(longitudes);
latitudes = latitudes(valid_indices);
longitudes = longitudes(valid_indices);
countries = countries(valid_indices);

tbl = tabulate(countries);
tbl = sortrows(tbl,1);
attack_counts = cell2mat(tbl(:, 2));
sorted_data = sortrows(tbl, -2);
% use only one pair of latitude and longitudes to represent a country
unique_countries = unique(countries); % Get unique country names

% Initialize arrays to store aggregated latitude and longitude data
avg_latitudes = zeros(size(unique_countries));
avg_longitudes = zeros(size(unique_countries));

% Loop through unique countries and calculate average latitude and longitude
for i = 1:numel(unique_countries)
    country = unique_countries{i};
    indices = strcmp(countries, country); % Find indices of data corresponding to the current country
    avg_latitudes(i) = mean(latitudes(indices)); % Calculate mean latitude
    avg_longitudes(i) = mean(longitudes(indices)); % Calculate mean longitude
end

% Create a new table or structure to store the unique country names along with their representative latitude and longitude
country_coordinates = table(unique_countries, avg_latitudes, avg_longitudes, 'VariableNames', {'Country', 'Latitude', 'Longitude'});

% Merge attack counts with country coordinates based on the country names
country_coordinates.attack_counts = attack_counts;
figure('Name',"Attack Location 2000-2021")
% Geographic plot
geobubble(country_coordinates.Latitude, country_coordinates.Longitude, country_coordinates.attack_counts);

% Customize the map properties

title('Attack Location 2000-2021');

%% Attack Location 1970-2000
data_early = data(data.iyear <= 2000, :);
% Extract latitude and longitude information
latitudes = data_early.latitude;
longitudes = data_early.longitude;
countries = data_early.country_txt;

% Remove missing or NaN values in latitude, longitude, and year
valid_indices = ~isnan(latitudes) & ~isnan(longitudes);
latitudes = latitudes(valid_indices);
longitudes = longitudes(valid_indices);
countries = countries(valid_indices);

tbl = tabulate(countries);
tbl = sortrows(tbl,1);
attack_counts = cell2mat(tbl(:, 2));
sorted_data = sortrows(tbl, -2);
% use only one pair of latitude and longitudes to represent a country
unique_countries = unique(countries); % Get unique country names

% Initialize arrays to store aggregated latitude and longitude data
avg_latitudes = zeros(size(unique_countries));
avg_longitudes = zeros(size(unique_countries));

% Loop through unique countries and calculate average latitude and longitude
for i = 1:numel(unique_countries)
    country = unique_countries{i};
    indices = strcmp(countries, country); % Find indices of data corresponding to the current country
    avg_latitudes(i) = mean(latitudes(indices)); % Calculate mean latitude
    avg_longitudes(i) = mean(longitudes(indices)); % Calculate mean longitude
end

% Create a new table or structure to store the unique country names along with their representative latitude and longitude
country_coordinates = table(unique_countries, avg_latitudes, avg_longitudes, 'VariableNames', {'Country', 'Latitude', 'Longitude'});

% Merge attack counts with country coordinates based on the country names
country_coordinates.attack_counts = attack_counts;

figure('Name',"Attack Location 1970-2000")
% Geographic plot
geobubble(country_coordinates.Latitude, country_coordinates.Longitude, country_coordinates.attack_counts);

% Customize the map properties

title('Attack Location 1970-2000');

%% Attack Per Continent 

% Consolidate continents as specified
data.region_txt(ismember(data.region_txt, {'Central Asia', 'East Asia', 'South Asia', 'Southeast Asia'})) = {'Asia'};
data.region_txt(ismember(data.region_txt, {'Middle East & North Africa', 'Sub-Saharan Africa'})) = {'Africa'};
data.region_txt(ismember(data.region_txt, {'Western Europe', 'Eastern Europe'})) = {'Europe'};
data.region_txt(ismember(data.region_txt, {'Central America & Caribbean'})) = {'North America'};
data.region_txt(ismember(data.region_txt, {'Australasia & Oceania'})) = {'Oceania'};

% Desired continents
desired_continents = {'Africa', 'Asia', 'Europe', 'North America', 'South America', 'Oceania', 'Antarctica'};

% Group data by year and continent, then count attacks
attack_counts = grpstats(data, ["iyear", "region_txt"]);

% Reshape the data for plotting
years = unique(attack_counts.iyear);
num_years = numel(years);
num_continents = numel(desired_continents);

attack_counts_reshaped = zeros(num_years, num_continents);
for i = 1:num_years
    for j = 1:num_continents
        idx = attack_counts.iyear == years(i) & strcmp(attack_counts.region_txt, desired_continents{j});
        if any(idx)
            attack_counts_reshaped(i, j) = attack_counts.GroupCount(idx);
        end
    end
end


figure('Name','Attack per Continent');

% Define x-axis positions for bars
x_positions = 1:num_years;

% Plot stacked bars using the specified x-axis positions
bar(x_positions, attack_counts_reshaped, 'stacked');

% Update x-axis ticks and labels to display years
xticks(x_positions);
xticklabels(string(years));
xlabel('Year');
ylabel('Number of Attacks');
title('Number of Attacks per Continent over the Years');
legend(desired_continents, 'Location', 'best');



%% Number of Deaths per Year
% Extract the year and number of deaths
years = data.iyear;
num_deaths = data.nkill;

% Remove missing or NaN values in the number of deaths and nonpositive values
valid_data = ~isnan(num_deaths) & num_deaths > 0 & ~isnan(years);

% Filter data based on valid entries
years = years(valid_data);
num_deaths = num_deaths(valid_data);

[counts, x] = hist(num_deaths, 10.^(0:0.1:3)); %Compute the density of degree intervals
figure
subplot(2,1,1)
[power] = logfit(x,counts,'loglog');
xlabel('Deaths');
ylabel('P(Deaths)');
title('Deaths distribution')
axis square
grid on


subplot(2,1,2)
plot(num_deaths)
axis tight
grid on
ylabel('Nunber of Deaths');
title('Deaths per attack')
xticks(1:5900:length(years))
xticklabels(years(1:5900:end))
xtickangle(30)


%% Deaths divided by Target types

% Extract relevant columns
years = data.iyear;
num_deaths = data.nkill;
target_types = data.targtype1_txt;
% Remove missing or NaN values in the number of deaths, nonpositive values, and invalid years
valid_data = ~isnan(num_deaths) & num_deaths > 0 & ~isnan(years);

% Filter data based on valid entries
years = years(valid_data);
num_deaths = num_deaths(valid_data);
target_types = target_types(valid_data);

% Get unique target types and years
unique_target_types = unique(target_types);
unique_years = unique(years);

% Initialize matrix to store number of deaths per target type per year
deaths_by_target_type = zeros(length(unique_years), length(unique_target_types));

% Calculate number of deaths per target type per year
for i = 1:length(unique_years)
    for j = 1:length(unique_target_types)
        deaths_by_target_type(i, j) = sum(num_deaths(years == unique_years(i) & strcmp(target_types, unique_target_types{j})));
    end
end

% Define a colormap with enough colors for unique_target_types
colors = hsv(length(unique_target_types)); % You can use any colormap you prefer

% Define x-axis positions for bars
x_positions = 1:num_years;

% Create stacked bar plot with different colors for each target type
figure('Name','Death per Target Type');
h = bar(x_positions, deaths_by_target_type, 'stacked');

% Assign colors to each set of bars representing a target type
for k = 1:length(unique_target_types)
    set(h(k),'FaceColor', colors(k,:));
end

xticks(x_positions);
xticklabels(string(unique_years));
xlabel('Year');
ylabel('Number of Deaths');
title('Number of Deaths by Target Type over the Years');
legend(unique_target_types, 'Location', 'bestoutside');

%% Number of Attacks per Target Type

% Extract relevant columns
target_types = data.targtype1_txt;
tab = tabulate(target_types);

% Extract target types and their corresponding counts
unique_target_types = tab(:, 1);
counts = cell2mat(tab(:, 2));

% Create a horizontal bar plot for the number of attacks per target type
figure('Name', 'Number of Attacks per Target Type');
bar(unique_target_types,counts);
xlabel('Number of Attacks');
ylabel('Target Type');
title('Number of Attacks per Target Type');
