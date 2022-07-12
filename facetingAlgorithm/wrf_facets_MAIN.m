%% Savanna Wolvin
% created: Aug. 19th, 2020
% edited: Jul. 11th, 2022

% SUMMARY
% Create facets for wrf gridpoints. You choose the sampling radius
% which determines the facet facing direction.

% INPUT
% wrf_file -    file path to the WRF files, using elevation/lat/lon data
% dist_grid -   Use this variable or dist_km to choose sampling radius. 
%               This variable chooses sampling radius by number of grid 
%               points around
% dist_km -     This variable chooses sampling radius by a kilometer
%               radius
% th -          From the calculated rise over run between each gridpoint, 
%               this is the maximum change in height between gridpoints to
%               be designated as flat
% edges -       The facing directions of the facets are binned, these are
%               the edges of the bins

% OUTPUT
% topography.png -      Plot of smoothed topography
% flat_terrain.png -    Plot of flat VS not flat terrain
% facet_??km.png -      Plot of facets with directions shaded
% wrf_facets_??.m -     .mat file with facets_orientation (datapoints 
%                       numbered by facing direction) and facets_labeled 
%                       (datapoints numbered by facet)



%% Clear workspace and add paths

%clear; close all;

% mapping functions
% https://www.eoas.ubc.ca/~rich/map.html
addpath('/uufs/chpc.utah.edu/common/home/strong-group4/savanna/m_map');

% export figure path
% https://www.mathworks.com/matlabcentral/fileexchange/23629-export_fig
addpath('/uufs/chpc.utah.edu/common/home/strong-group4/savanna/export');

% countries shapefiles
% http://www.naturalearthdata.com/
addpath('/uufs/chpc.utah.edu/common/home/strong-group4/savanna/countries_shapefile');




%% Variable presets

% wrf file; where the terrain data will come from
fdir    = "/uufs/chpc.utah.edu/common/home/strong-group6/HIMAT/WRF_output/2000/output/";
fname   = "wrfout_d03_2000-01-01_00:00:00"; 
wrf_file    = fdir + fname;

dist_grid   = 4; % distance from location by gridpoints
dist_km     = 25; % distance from location by kilometers
th      = 5; % max angle for terrain to be flat. (degrees)
edges   = -180:45:180; % defining bin sizes, this one had 8 bins

% directory to save the figures and the .mat file holding the facets
save_dir = '/uufs/chpc.utah.edu/common/home/strong-group4/savanna/wrf_facet_data/';

% set TRUE if you want grid point and test/rule outcome displayed
outcome = true;


%% Pull WRF data

disp("Load terrain data fro WRF model...");

elev = ncread(wrf_file, 'HGT', [1 1 1], [inf inf 1]);
lati = ncread(wrf_file,'XLAT',[1 1 1],[inf inf 1]); lati = double(lati); 
loni = ncread(wrf_file,'XLONG',[1 1 1],[inf inf 1]); loni = double(loni); 

latlim = [min(lati(:)) max(lati(:))];
lonlim = [min(loni(:)) max(loni(:))]; 




%% Elevation Calculations

% smooth elevation by guassian, w/ STD of sigma
elev = imgaussfilt(elev,40/3/4);  % 3-sigma is 30 km, so sigma = 30/3/4km

% calculate direction topo is facing
[dx, dy] = gradient(elev);
dy = -dy;
dx = -dx; dy = -dy; % flip orientation

% angle of the direction the topo faces. (Degrees)
angs = atan2(dy(:),dx(:))*180/pi;
angs = reshape(angs,size(lati)); % back to lat-lon grid

% determine what gridpoints are flat. FLAT EARTH
flat = abs(dx)<th & abs(dy)<th; % slope<=prctile(slope,60);

% creates histogram of angles within each bin
[~,comp] = histc(angs,edges); % comp = matrix of bin num for each gridpoint

% mesh grid for defining distances from gridpoint
[xx,yy] = meshgrid(1:size(elev,2),1:size(elev,1));
close all;  




%% Create a map of the topography

disp("Plotting the smoothed topography to be used...");


fig1 = figure; %('visible', 'off');
m_proj('lambert', 'lon', lonlim, 'lat', latlim);
m_contourf(loni, lati, elev, 'edgecolor', 'none');
m_coast();
colormap([m_colmap('green') ;m_colmap('bland')]);
M=m_shaperead('ne_50m_admin_0_countries'); 
for k=1:length(M.ncst)
     m_line(M.ncst{k}(:,1),M.ncst{k}(:,2), 'color', [0 0 0]); 
end
title('Topography of High Mountain Asia');
m_grid('tickdir', 'out', 'yaxislocation', 'right', 'xaxislocation', 'bottom');
colorbar('southoutside');

fig_name = save_dir + "topography.png";
export_fig(fig1, fig_name{:}, '-transparent', '-r300', '-nocrop')




%% Map of Flat VS Not Flat

disp("Plotting where the terrain is flat vs not flat...");

fig2 = figure;
m_proj('lambert', 'lon', lonlim, 'lat', latlim);
m_coast('patch', [.7 .7 .7]);
m_contourf(loni, lati, real(flat), 'edgecolor', 'none');
M=m_shaperead('ne_50m_admin_0_countries'); 
for k=1:length(M.ncst)
     m_line(M.ncst{k}(:,1),M.ncst{k}(:,2), 'color', [0 0 0]); 
end
cmap=[0, 0.6, 0.011 ; 1, 1, 1];
colormap(cmap);
title('Mountainous Regions');
m_grid('tickdir', 'out', 'yaxislocation', 'right', 'xaxislocation', 'bottom');
colorbar('southoutside',  'Ticks', [0 1], 'tickLabels',{'Not Flat','Flat'});

fig_name = save_dir + "flat_terrain.png";
export_fig(fig2, fig_name{:}, '-transparent', '-r300', '-nocrop')




%% Defining the Facets for Each Gridpoint

disp("Determining grid point orientations...");

facets = NaN(size(elev)); % empty facets matrix

for i = 1:size(comp,1)      % lon loop
    for j = 1:size(comp,2)  % lat loop (PARFOR)
        if outcome; disp("[" + string(i) + ", " + string(j) + "]"); end
        
        % FIRST: check if facet is flat. If flat, skip loop
        if flat(i,j) || isnan(elev(i,j))
            facets(i,j) = 9;
            continue;
        end
        
        
        % NOT FLAT so pull locations around that gridpoint

        % Distance by gridpoint
%         use = sqrt((xx-xx(i,j)).^2 + (yy-yy(i,j)).^2)<dist_grid; % locations around gridpoint

        % Distance by arclength
        [arclen, ~] = distance(lati(i,j), loni(i,j), lati, loni, 'degrees');
        use = deg2km(arclen) < dist_km;
        
        % bin the locations around gridpoint [num of elements in bin, bin num]
        [c,bin] = histc(comp(use),1:8); 
        fu = flat(use);     % locations around gridpoint that are flat
        ac = [c (1:8)'];    % rank and apply rule. 
        sumc = sum(c);      % num of datapoints within 4 indecies
        ac = flipud(sortrows(ac,1)); % decreasing sort, most to least
        ac(:,1) = ac(:,1)/sumc; % precent of each bin found within 4 indecies
        fracs = ac(:,1);    % fraction
        comps = ac(:,2);    % bin number
        flb = [];           % will contain percent flat  
         
        facets = gibson_tree(i, j, fracs, comps, facets, flb, bin, fu, ac, outcome);
        
    end
end    


%% Mapping the facets
% a place to re-load the facets in case you dont want to re-run the script
% load("/uufs/chpc.utah.edu/common/home/strong-group4/savanna/wrf_facet_data/wrf_facets_25");
%%

disp("Plotting orientations of the grid points...");

fig3 = figure('Position', [10 10 800 600]);%('visible', 'off');

m_proj('sinusoidal', 'lon', lonlim, 'lat', latlim);
hold on 
m_coast('patch', [1 1 1], 'edgecolor', 'none');
m_contourf(loni, lati, facets, 'edgecolor', 'none');
M=m_shaperead('ne_50m_admin_0_countries'); 
% for k=1:length(M.ncst)
%      m_line(M.ncst{k}(:,1),M.ncst{k}(:,2), 'color', [105, 68, 28]./255,'linewi',1.5); 
% end

for k=1:length(M.ncst)
     m_line(M.ncst{k}(:,1),M.ncst{k}(:,2), 'color', [105, 68, 28]./255,'linewi',2); 
end

M=m_shaperead('ne_50m_coastline'); 
for k=1:length(M.ncst)
     m_line(M.ncst{k}(:,1),M.ncst{k}(:,2), 'color', [0 73 110]./255,'linewi',2.5); 
end


cmap = rot90(0:36.42857:255)/255;
cmap = flipud([cmap cmap cmap]);
cmap = [cmap; 1 0.92549 0.76078];

colormap(cmap);
title('Facets of High Mountain Asia','fontsize', 16);
m_grid('tickdir', 'in', 'linest','none', 'yaxislocation', 'right', 'xaxislocation', 'bottom','fontsize', 12);
cbr = colorbar('southoutside', 'Ticks', [1.4 2.4 3.2 4.1 5 5.9 6.8 7.7 8.6], 'tickLabels',...
    {'SSE','ESE','ENE','NNE','NNW','WNW','WSW','SSW','FLAT'},'fontsize', 12, 'TickLength', 0.0001);

fig_name = save_dir + "facet_" + dist_km + "km.png";
export_fig(fig3, fig_name{:}, '-transparent', '-r300', '-nocrop')



%% this loop adds unique number labels to each facet

disp("Labeling the facets...");

max_label = 0;
bw = zeros(size(facets));
for i=1:9 % changed to 9 from 8 to get flat facets
    % labels each facet by number from a certain direction, like a paint by number :)
    directional_facet = bwlabel(facets==i); 
    % logical matrix identifying locations of facets in that direction
    found_facets = directional_facet>0;
    % set bw with previous bw, current directional_facet, and offset with st
    bw(found_facets) = bw(found_facets) + max_label + directional_facet(found_facets); 
    max_label = max(bw(:));         % reset st to the highest current facet label
end
facets_labeled = bw;
num_of_facets = unique(facets_labeled(:));




%% save the facet file with the array holding the orienation of each gridpoint
% and the array holding the numbered facets

disp("Saving facets in MAT file...");

facets_orientation = facets;
save(save_dir + "wrf_facets_" + dist_km,'facets_orientation','facets_labeled'); 


