function [surveymap, bins, shortest] = plotSurvey(bin_width, lat_dist, lon_dist, track, waypoint, image_file_p, world_file_p)
%PLOTSURVEY creates a matrix using the gps track which was imported using
% import_File.m and computed using surveyDim.m. PLOTSURVEY scores each grid
% cell based on whether the cell was visited, viewed or missed. A visited
% cell gets given a value of 1, a viewed cell 0.5 and a missed cell 0.
% Based on the cell sizes, this gives an indication of how well an area was
% surveyed with the aim being to have no missed cells. The figure is 
% plotted using the function PCOLOR and has axes values in meters according 
% to the UTM coordinate system. Observations can also be included, this
% will give a grid cell the value of 2 and this is based on the waypoint
% input variable.
%
% Input arguments
% ----------------
% bin_width = size of each grid cell
%             A bin_width of 2 creates 2x2 m cells (4 m^{2}).
% lat_dist = latitudinal distance obtained from surveyDim.m
% lon_dist = longitudinal distance obtained from surveyDim.m
%            Both lat_dist and lon_dist, together with bin_width, are used
%            to set the size of the cells. 
%            i.e. lon_dist = 397; bin_width = 2; 
%                 lon_bin = lon_dist/bin_width;
%                 lon_bin = 199;
%                 Therefore there are 199, 2 m bins along the x-axis. 
% track = m x 2 [lon, lat] array
%         Fills the grid cells with data. Each lat/lon pair which fall in a
%         grid cell give that cell a count. If a cell has at least one
%         lat/lon pair, that cell has been 'visited' by a surveyor.
%
% waypoint = set of lat/lon coordinates. 
%            Is added to data generated from track and creates a map with
%            observations.
%
% image file, world file = A map image (.jpg or .png) (from Google Earth) 
%                          and it's associated world file (.jgw or .pgw). 
%                          These are used my the picplot.m function which 
%                          plots the survey track grid over the map. Please
%                          see documentation on how to generate a world
%                          file for a google earth image.
%
% Output variables
% -----------------
% surveymap = gridded dataset which represents visited, viewed and missed
%             grid cells with 1, 0.5 and 0 resepctivley.
% bins = [lon_bins, lat_bins]
%        Indicates the size of the matrix and the scaling of the bins, also
%        used to calculate efficacy coefficients with surveyEfficacy.m.
% shortest = shortest distance from observation to missed cell.
%
% Example 1: Small or overgrown survey area - bin_width = 1
% --------
% [surveymapcsv, binscsv, shortest] = plotSurvey(1, lat_dist_csv, lon_dist_csv, track_csv);
% 
% Example 2: large or open survey area - bin_width = 5
% --------
% [surveymapcsv, binscsv, shortest] = plotSurvey(5, lat_dist_csv, lon_dist_csv, track_csv);
%

% Example 3: Adding a map image. For the Bien Donne site.
% [surveymap, bins] = plotSurvey(1, lat_dist, lon_dist, t, 'BDF.jpg', 'BDF.jgw');

%% Google Earth photos

if ~exist('track', 'var') || ~exist('waypoint', 'var')
    error('goodnessOfSurvey:plotSurvey:missingVariable','Input variable track/waypoint is missing')
end

if ~exist('bin_width', 'var') || ~exist('lat_dist', 'var') || ~exist('lon_dist', 'var')
    error('goodnessOfSurvey:plotSurvey:missingVariable','Input variable bin_width/lat_dist/lon_dist is missing')
end 

if ~exist('image_file_p', 'var')
    image_file_p = {};
end
if ~exist('world_file_p', 'var')
    world_file_p = {};
end

%% Set bin size

lon_bin = double(lon_dist/bin_width);
lat_bin = double(lat_dist/bin_width);
bins = [lon_bin, lat_bin];

waypoint = [waypoint(:,2), waypoint(:,1)];

% Create gridded GPS track: indicates which Xm^2 areas are covered by
% a surveyor. Data stored as a matrix.
figure('Visible', 'off');
[H, bc] = hist3(track, bins);
hist3(track, bins)
caxis([0, 1])
cmap = jet(256);
cmap(1,:) = 1;
colormap(cmap)
colorbar 
xlabel('Longitude'); 
ylabel('Latitude');
set(gcf,'renderer','opengl');
set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
view(2)
%Waypoints
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Visible', 'off');
[Hw] = hist3(waypoint, 'Ctrs', bc);
hist3(waypoint, bc)
caxis([0, 5])
cmap = jet(256);
cmap(1,:) = 1;
colormap(cmap)
colorbar 
xlabel('Longitude'); 
ylabel('Latitude');
set(gcf,'renderer','opengl');
set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
view(2)
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make all values greater than 1 equal 1
H(H >= 1) = 1;

for xx = 4:(lon_bin-3)
    for yy = 4:(lat_bin-3)
        if H(xx,yy) >= 1
            if H(xx+1,yy) == 0
                H(xx+1,yy) = 0.5;
            end
            if H(xx+2,yy) == 0
                H(xx+2,yy) = 0.5;
            end
            if H(xx+3,yy) == 0
                H(xx+3,yy) = 0.5;
            end
            if H(xx-1,yy) == 0
                H(xx-1,yy) = 0.5;
            end
            if H(xx-2,yy) == 0
                H(xx-2,yy) = 0.5;
            end
            if H(xx-3,yy) == 0
                H(xx-3,yy) = 0.5;
            end            
            if H(xx,yy+1) == 0
                H(xx,yy+1) = 0.5;
            end
            if H(xx,yy+2) == 0
                H(xx,yy+2) = 0.5;
            end
            if H(xx,yy+3) == 0
                H(xx,yy+3) = 0.5;
            end
            if H(xx,yy-1) == 0
                H(xx,yy-1) = 0.5;
            end
            if H(xx,yy-2) == 0
                H(xx,yy-2) = 0.5;
            end
            if H(xx,yy-3) == 0
                H(xx,yy-3) = 0.5;
            end
        end
    end
end

H_obs = H + Hw;

% Make 0 values outside of survey area NaN (L to R: essentially top to bottom as before transpose)
for xa = 1:lon_bin
    for xb = 1:lat_bin
        if H(xa,xb) == 0
            H(xa,xb) = NaN;
            H_obs(xa,xb) = NaN;
        elseif H(xa,xb) > 0
           break
        end  
    end
end
% Make 0 values outside of survey area NaN (R to L: essentially bottom to top)
for ya = lon_bin:-1:1
    for yb = lat_bin:-1:1
        if H(ya,yb) == 0
            H(ya,yb) = NaN;
            H_obs(ya,yb) = NaN;
        elseif H(ya,yb) > 0
           break
        end  
    end
end

% Transpose so orientation of matrix is correct 
surveymap = H'; 
obsmap = H_obs';
obsmap(obsmap >= 2) = 2;
% Make 0 values outside of survey area NaN (L to R)
for xc = 1:lat_bin
    for xd = 1:lon_bin
        if surveymap(xc,xd) == 0
            surveymap(xc,xd) = NaN;
            obsmap(xc,xd) = NaN;
        elseif surveymap(xc,xd) > 0
           break
        end  
    end
end
% Make 0 values outside of survey area NaN (R to L)
for yc = lat_bin:-1:1
    for yd = lon_bin:-1:1
        if surveymap(yc,yd) == 0
            surveymap(yc,yd) = NaN;
            obsmap(yc,yd) = NaN;
        elseif surveymap(yc,yd) > 0
           break
        end  
    end
end
            
% Create X and Y axis: latitude and longitude need to be the same length
% as the dimensions of matrix H1_corrected.
%lat_ax = lat_min:((lat_max - lat_min)/(lat_bin-1)):lat_max;
%lon_ax = lon_min:((lon_max - lon_min)/(lon_bin-1)):lon_max;

% Axes in UTM

tracks_utm = ll2utm(track(:,2), track(:,1));
lat_utm_max = max(tracks_utm(:,2));
lat_utm_min = min(tracks_utm(:,2));
lon_utm_max = max(tracks_utm(:,1));
lon_utm_min = min(tracks_utm(:,1));
lat_ax = lat_utm_min:((lat_utm_max - lat_utm_min)/(lat_bin-1)):lat_utm_max;
lon_ax = lon_utm_min:((lon_utm_max - lon_utm_min)/(lon_bin-1)):lon_utm_max;

% Axes in degrees

%lat_ax_deg = min(track(:,2)):((max(track(:,2))-min(track(:,2)))/(lat_bin-1)):max(track(:,2));
%lon_ax_deg = min(track(:,1)):((max(track(:,1))-min(track(:,1)))/(lon_bin-1)):max(track(:,1));


% Plot matrix as a surface plot. Shows grids where a surveyor walked 
% as equivalent to 1, and neighbouring grids which were in viewing range
% as 0.5. Observations where grid value greater than 1.
figure('Visible', 'On');
if isempty(image_file_p)
    disp("No image file selected.")
elseif ~regexp(image_file_p, regexptranslate('wildcard', '*.jpg')) 
    OR ~regexp(image_file_p, regexptranslate('wildcard', '*.png'))
    error('goodnessOfSurvey:plotSurvey:fileExtensionError','image_file_p must be .jpg or .png')
elseif ~regexp(world_file_p, regexptranslate('wildcard', '*.jgw')) 
    OR ~regexp(world_file_p, regexptranslate('wildcard', '*.pgw'))
        error('goodnessOfSurvey:plotSurvey:fileExtensionError','world_file_p must be .jgw or .pgw')
else 
    picplot(image_file_p, world_file_p);
    hold on
end
pc = pcolor(lon_ax, lat_ax, obsmap);
caxis([0, 2])
cmap = bluewhitered(4);
cmap(1,:) = [1 0 0];
cmap(2,:) = [1 0.5 0];
cmap(3,:) = [0 0 1];
cmap(4,:) = [0 1 0];
colormap(cmap)
cb = colorbar;
cb.FontSize = 16;
cb.YTick = [2/8 6/8 10/8 14/8];
cb.YTickLabel = {'Missed', 'Viewed', 'Visited', 'Observation'};
xlabel('UTM_{x} [m]', 'FontSize', 16); 
ylabel('UTM_{y} [m]', 'FontSize', 16);
hAx = gca;
hAx.XAxis.FontSize = 14;
hAx.YAxis.FontSize = 14;
hAx.YAxis.Exponent=0;
hAx.XAxis.Exponent=0;
ytickformat('%d')
xtickformat('%d')
set(hAx,'XMinorTick','on','YMinorTick','on')
set(pc, 'EdgeColor', 'none');
%xlim([min(lon_ax)-5, max(lon_ax)+5])
%ylim([min(lat_ax)-5, max(lat_ax)+5])
str = strcat("Grid cell area = ", num2str(bin_width), 'm^{2}');
annotation('textbox',[0.751 0.495 0.5 0.5],'String',str,'FitBoxToText','on', 'FontSize', 14);
hold off

% Get shortest distances from observations to missed
[obs_row, obs_col] = find(obsmap == 2);
[mis_row, mis_col] = find(obsmap == 0);

shortest = zeros(length(obs_row),1);
permissedpoint = zeros(length(mis_row),1);
for d = 1:length(obs_row)
    for dd = 1:length(mis_row)
        permissedpoint(dd) = sqrt((obs_col(d) - mis_col(dd)).^2 + (obs_row(d) - mis_row(dd)).^2);
    end
    shortest(d) = min(permissedpoint);
    permissedpoint(:) = 0;
end

shortest = min(shortest);
disp(shortest)
        
end
