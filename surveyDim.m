function [track, lat_dist, lon_dist] = surveyDim(varargin)
%SURVEYDIM Calculates the latitudinal and longitudinal distances of a 
% survey area given the latitudes and longitudes as input arguments. This 
% function is modeled on the haversine equation which uses latitude and 
% longitude (as radians) to calculate distances in meters. The purpose of 
% this calculation is so the bin widths can be correctly sized for plotting. 
% If multiple files are imported using import_File, this function will also
% concatenate them into a single m x 2 [lon, lat] array. 
%
% Input arguments
% ----------------
% varargin = any number of m x 2 [lat, lon] arrays
%            Most surveys will consist of multiple files, list all imported
%            [lat, lon] arrays which were imported using 'import_File'.
%
% Output variables
% -----------------
% track = m x 2 [lon, lat] array
% lat_dist = latitudinal distance
% lon_dist = longitudinal distance
%
% Example 1
% ----------
% [track_csv, lat_dist_csv, lon_dist_csv] = surveyDim(trackcsv);
% 'track_csv' [lon, lat] and 'trackcsv' [lat, lon] will be the same size.
%
% Example 2
% ----------
% [track_gpx, lat_dist_gpx, lon_dist_gpx] = surveyDim(trackgpx1, trackgpx2, ..., trackgpxN);
% 'track_gpx' [lon, lat] will be size of the sum of the lengths of 
% 'trackgpx1', 'trackgpx2', ..., 'trackgpxN' [lat, lon].
                
%% Compile all imported tracks into one m x 2 [lon, lat] array
%nargin
if nargin < 1
    error("goodnessOfSurvey:surveyDim:notEnoughInputArguments", "Atleast one set of lat/lon coordinates are reauired.")
end

if nargin == 1
    track = [varargin{1}(:,2), varargin{1}(:,1)];
elseif nargin > 1
    arrlen = zeros(1, nargin);
    for x = 1:nargin    
        arrlen(x) = length(varargin{x});
    end
    track = zeros(sum(arrlen), 2);
    for n = 1:nargin
        track(1:length(varargin{n}), :) = [varargin{n}(:,2), varargin{n}(:,1)];
    end
end

%% Haversine Calculation
% Convert to radians
lat_rad = track(:,2).*(pi/180);
lon_rad = track(:,1).*(pi/180);

% Delta lat and delta lon
dlon_rad = max(lon_rad) - min(lon_rad);
dlat_rad = max(lat_rad) - min(lat_rad);

% Haversine Formula
% Latitudinal distance
Alat = sin(dlat_rad/2).^2 + cos(max(lat_rad)) * cos(min(lat_rad)) * sin(0/2).^2;
Clat = 2 * atan2(sqrt(Alat), sqrt(1-Alat));
Dlat = (6371 * Clat) * 1000;
lat_dist = int16(Dlat);

% Longitudinal distance
Alon = sin(0/2).^2 + cos(max(lat_rad)) * cos(min(lat_rad)) * sin(dlon_rad/2).^2;
Clon = 2 * atan2(sqrt(Alon), sqrt(1-Alon));
Dlon = (6371 * Clon) * 1000;
lon_dist = int16(Dlon);

end

    