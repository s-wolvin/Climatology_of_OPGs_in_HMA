%% Savanna Wolvin
% created: Mar. 11th, 2021
% edited: Nov. 24th, 2021

% SUMMARY

% INPUT
% wrf_dir       - Directory holding the WRF output files
% save_dir      - Directory to save the precip .mat files
% yr_initail    - Inital year of the WRF output files
% yr_final      - Final year of the WRF output files
% hours         - All hours to check

% OUTPUT
% plot_prcntOfTemp_YR_INITIAL_YR_FINAL.png - plot of minimum and maximum
% temperature frequencies over a 24 hour period of all mountainous
% gridpoints


%% clear workspace
clear; close all;

addpath('/uufs/chpc.utah.edu/common/home/u1324060/exp_fig');




%% variable presets

% WRF output data location
wrf_dir = "/uufs/chpc.utah.edu/common/home/strong-group6/HIMAT/WRF_output/";

% Save folder
save_dir = "/uufs/chpc.utah.edu/common/home/strong-group4/savanna/precip_wrf/";

% start and end year
yr_initial = 2001;
yr_final = 2015;

yr = yr_initial:yr_final;

hours = 0:23;

num_min = zeros(length(hours),1);
num_max = zeros(length(hours),1);




%% map presets
disp('Load map data...');

% wrf file
fdir    = "/uufs/chpc.utah.edu/common/home/strong-group6/HIMAT/WRF_output/2000/output/";
fname   = "wrfout_d03_2000-01-01_00:00:00"; 
wrf_file    = fdir + fname;

lati = ncread(wrf_file,'XLAT',[1 1 1],[inf inf 1]); lati = double(lati); 
loni = ncread(wrf_file,'XLONG',[1 1 1],[inf inf 1]); loni = double(loni); 
elev = ncread(wrf_file, 'HGT', [1 1 1], [inf inf 1]);
elev = elev';


%% Elevation Calculations
disp('Determine Flat Facets...');
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
th   = 5;
flat = abs(dx)<th & abs(dy)<th; % slope<=prctile(slope,60);




%% Loop through all years of data and pick hours in which minimum and 
%   maximum temperatures

disp('Calculate Histogram of Frequencies...');
for yrx = yr
    file_path = wrf_dir + num2str(yrx) + "/output/";
    di = dir(file_path + "wrfout_d03*"); % pull names of all files    

    temp = [];
    wtm = [];

    % loop through each file of that year and create variable of all temps
    % over one year
    for i=1:length(di)
        % disp file name
        fname = di(i).name;
        disp(fname);

        % read the time values for one month of the year & convert to nums
        wtmx = ncread(file_path + fname,'Times')';
        wtm = cat(1,wtm,datenum(wtmx));
        
        % pull surface temperatures
        tempx = ncread(file_path + fname,'T2');
        % *** LON X LAT X TIME ***
        
        % *** TIME X LAT X LON ***
        tempx = permute(tempx,[3 2 1]);
        temp = cat(1, temp, tempx);
        clear tempx
        
    end
    
    tindex = yrx == year(wtm);
    temp = temp(tindex,:,:); % pull only temp from days of that year
    temp(:,flat == 1) = NaN; % set temp at flat locations to NaN
    wtm = wtm(tindex);
    days = datestr(wtm, 'mm/dd/yyyy'); % create list of date strings
    daysLoop = unique(days, 'rows'); % pull all unique days (filter out hours)

    % loop through each day to count number of times a min/max temperature
    % occurs at each non-flat gridpoint
    for dayx = 1:length(daysLoop)
        % pull temps of single day & calculate index of 
        dayx_index = datenum(daysLoop(dayx, :), 'mm/dd/yyyy') == datenum(days(:,:), 'mm/dd/yyyy');
        [~, minIdx] = min(temp(dayx_index,:,:));
        [~, maxIdx] = max(temp(dayx_index,:,:));
        
        minIdx = squeeze(minIdx); % squeeze out extra dimension
        maxIdx = squeeze(maxIdx);
        
        minIdx(flat == 1) = NaN; % remove flat facets
        maxIdx(flat == 1) = NaN;
        
        minIdx = sum(hist(minIdx, 1:24),2); % count number of gridpoints
        maxIdx = sum(hist(maxIdx, 1:24),2); % with min/max temps occurring
        
        num_min = sum([num_min, minIdx], 2); % add histograms together
        num_max = sum([num_max, maxIdx], 2);
    end
     
    
    clear temp
    clear wtm
end




%% Plot frequncies
disp('Plot...');

fig2 = figure();

totalmin = sum(num_min,'all');
prcntmin = num_min ./ totalmin;

totalmax = sum(num_max,'all');
prcntmax = num_max ./ totalmax;

p1 = plot(0:23,prcntmin(:,:));
p1.LineWidth = 2;
hold on;
p2 = plot(0:23,prcntmax(:,:));
p2.LineWidth = 2;
xticks([0:2:23]);
xticklabels(["00","02","04","06","08","10","12","14","16","18","20","22"]);
xlim([0,23]);
ylabel("Frequency (%)",'Color', 'k');
yticks([0:0.10:0.60]);
yticklabels(["0","10","20","30","40","50","60","70","80"]);
xlabel("Hour (UTC)");
ylim([0,0.50]);

grid on
title("Frequency of Temperature Event at Each Hour");
legend(["Temp. Min.","Temp. Max."]);

fig_name = save_dir + "plot_prcntOfTemp_" + string(yr_initial) + "_" + string(yr_final) + ".png";
export_fig(fig2, fig_name{:}, '-transparent', '-r300')



