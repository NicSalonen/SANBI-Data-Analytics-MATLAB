function [track, filepath] = importSurvey(filepath, filename, numberOfLogs, Filetype)
%IMPORTSURVEY Imports lat/lon data from a CSV/GPX file of GPS track data. The user
% specifies the directory FILEPATH and file FILENAME as well as being able
% to specifiy the number of track logs NUMBEROFLOGS which need to be 
% imported if the file is a .gpx file. It is important to know the details
% of the file prior to importing, these can be found using mapping software
% such as the open source tool Garmin BaseCamp. 
%
% Input arguments
% ----------------
% filepath = path to working directory
%            Working directory should contain functions and data unless 
%            function filepaths' are added manually. 
% filename = name of CSV/GPX file
%            Must include the file extension (.csv or .gpx). 
% numberOfLogs = number of track logs/waypoints in a .gpx file
%                Imports the specified range of entries which can be 
%                stored in a .gpx file. To import all track logs/waypoints,
%                numberOfLogs must cover the entire index range. The upper
%                limit can exceed the range. The default = 1 and will only 
%                import the first track log or first waypoint.
% 
% Filetype = 'track' (default) or 'waypoint'
%            Specifies if .gpx file is a waypoint file. Default file type 
%            is track.
%
% Output variables
% -----------------
% track =  m x 2 [lat, lon] array
% filepath = the input filepath
%
% Examples
% ---------
% Output track and filepath
% --------------------------
% [trackcsv, pathcsv] = importSurvey("F:\PES2020", "Track Acacia cultriformis and fimbriata Makhanda Grey Dam 20201026-28.csv");
% [trackgpx, pathgpx] = importSurvey("F:\PES2020", "Track Acacia cultriformis and fimbriata Makhanda Grey Dam 20201026-28.gpx", 1:9);
% [trackwpt, pathwpt] = importSurvey("F:\PES2020", "Waypoint Acacia cultriformis and fimbriata Makhanda Grey Dam 20201026-28.gpx", 1:500, 'waypoint');
%
% Output track only
% ------------------
% trackcsv = importSurvey("F:\PES2020", "Track Acacia cultriformis and fimbriata Makhanda Grey Dam 20201026-28.csv");
% trackgpx = importSurvey("F:\PES2020", "Track Acacia cultriformis and fimbriata Makhanda Grey Dam 20201026-28.gpx", 1:9);
% trackwpt = importSurvey("F:\PES2020", "Waypoint Acacia cultriformis and fimbriata Makhanda Grey Dam 20201026-28.gpx", 1:500, 'waypoint');
%% Addpath to folder and subfolders
filepath = char(filepath);
addpath(genpath(filepath));
cd(filepath);

%% Import lat/lon from CSV/GPX
% Defaults
if ~exist('Filetype', 'var')
    Filetype = 'track';
    warning('goodnessOfSurvey:importSurvey:variableDefaultUsed', 'Filetype defaulted to track, if this is a waypoint file, set Filetype to waypoint.') 
end

if ~exist('numberOfLogs', 'var') && contains(filename, '.gpx')
    numberOfLogs = 1;
    warning('goodnessOfSurvey:importSurvey:variableDefaultUsed', 'GPX files may contain multiple track logs. numberOfLogs defaulted to 1, some track logs from input variable filename may be missing. numberOfLogs can be set from 1:n to obtain all track logs.') 
end
%
filename = char(filename);

if contains(filename, '.csv')
    trackcsv = readtable(filename);
    track = [trackcsv.lat, trackcsv.lon];
elseif contains(filename, '.gpx')
    trackgpx = gpxread(filename, 'FeatureType', Filetype, 'Index', numberOfLogs);
    track = [trackgpx.Latitude; trackgpx.Longitude];
    track = track';
    track = track(all(~isnan(track),2),:);

else
    error('goodnessOfSurvey:importSurvey:fileExtensionError', 'importSurvey only accepts .csv or .gpx file formats')
end
end