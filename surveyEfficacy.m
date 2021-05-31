function [efficacy] = surveyEfficacy(track, surveymap, bins, shortestDist)
%SURVEYEFFICACY Computes the efficacy of the survey.
% The efficacy coefficients are based on the areas which are visited,
% viewed and missed. The total area which is considered to be surveyed is
% the sum of the visited and viewed area. The percentage of area which is
% surveyed is computed using the bin_area which is then multiplied by the
% count of cells for each unique value (1, 0.5 and 0).
%
% Input arguments
% ----------------
% track = m x 2 [lon, lat] array
%         Used to obtain the latidudinal and longitudinal distances.
% surveymap = data grid
%             Contains information regarding surveyed status of each cell.
% bins = [lon_bin, lat_bin]
%        The number of longitudinal and latitudinal bins. The latitudinal
%        and longitudinal distances are divided by lon_bin and lat_bin
%        respectively to give an accurate estimation of a single cell's area.
%
% shortestDist = variable output from plotSurvey
%                Is added to summary table.
% Output variables
% -----------------
% efficacy = 1 x 4 table 
%            Contains information about the percentage of area which were
%            surveyed, visited, viewed and missed. The percentage of area
%            which was surveyed is the sum of the areas which were visited
%            and viewed, as well as the shortest distance from a plant
%            observation to a missed grid cell.
%
% Example
% --------
% efficacycsv = surveyEfficacy(track_csv, surveymapcsv, binscsv, shortDist);
% Output (comes with headings): 66.2241806664831 | 44.3018452217020 | 21.9223354447811 |
% 33.7758193335169 | 5.1

%% Generate Survey Efficacy coefficient
% Area of each grid cell/bin = bin_width^2
% Calculate precise bin area

utm_co = ll2utm(track(:,2), track(:,1));

bin_height = ((max(utm_co(:,2)) - (min(utm_co(:,2))))/bins(2));
bin_length = ((max(utm_co(:,1)) - (min(utm_co(:,1))))/bins(1));
bin_area = bin_length*bin_height;

area_missed = length(surveymap(surveymap == 0)).*bin_area;
area_viewed = length(surveymap(surveymap == 0.5)).*bin_area;
area_visited = length(surveymap(surveymap == 1)).*bin_area;
area_total = area_missed + area_viewed + area_visited;

% Percentage of area surveyed/viewed/missed
perc_missed = (area_missed/area_total)*100;
perc_visited = (area_visited/area_total)*100;
perc_viewed = (area_viewed/area_total)*100;
perc_surveyed = ((area_viewed+area_visited)/area_total)*100;

efficacy = table(perc_surveyed, perc_viewed, perc_visited, perc_missed, shortestDist, 'VariableNames', {'Surveyed', 'Viewed', 'Visited', 'Missed', 'Shortestdistmeters'});
end
