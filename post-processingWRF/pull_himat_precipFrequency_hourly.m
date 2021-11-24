%% Savanna Wolvin
% created: Feb. 20th, 2021
% edited: Nov. 24th, 2021

% SUMMARY

% INPUT
% wrf_dir -     directory holding the WRF output files
% save_dir -    directory to save the precip .mat files
% yr_initail -  inital year of the WRF output files
% yr_final -    final year of the WRF output files

% OUTPUT


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
yr_final = 2001;

yr = yr_initial:yr_final;

hours = 0:23;

precipSum = NaN(length(hours), length(yr));

frequency = zeros(length(hours),3);




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




%% Loop through all years of data and pick hours in which differing 
% precipiation frequencies occurrs

for yrx = yr
    file_path = wrf_dir + num2str(yrx) + "/output/";
    di = dir(file_path + "wrfout_d03*"); % pull names of all files    

    p = [];
    wtm = [];

    % loop through each file of that year and create variable of all precip
    % over one year
    for i=1:length(di)
        % disp file name
        fname = di(i).name;
        disp(fname);

        % read the time values for one month of the year & convert to nums
        wtmx = ncread(file_path + fname,'Times')';
        wtm = cat(1,wtm,datenum(wtmx));
        
        % pull Accumulated total grid scale precipitation, Microphysics
        % scheme
        pr1 = ncread(file_path + fname,'RAINNC');
        % pull accumulated total cumulus precipitation, Cumulus scheme
        pr2 = ncread(file_path + fname,'RAINC');
        % add precip together *** LON X LAT X TIME ***
        pr = pr1+pr2;
        
        % *** TIME X LAT X LON ***
        pr = permute(pr,[3 2 1]);
        p = cat(1, p, pr);
        clear pr
        
    end

    % turning the cumulative values of the precip and turn it into hourly values
    p(1:end-1,:,:) = p(2:end,:,:) - p(1:end-1,:,:);
    p = p(1:end-1,:,:);
    p(:,flat == 1) = 0;
    wtm = wtm(1:end-1);
   
    % loop through each day to count number of times < 1 mm/hr, < 3mm/hr
    % and >= 1 mm/hr, and > 3 mm/hr occurs at each non-flat gridpoint
    for hoursx = hours
        % find hours wanted
        timex = hour(wtm)==hoursx; 
        
        frequency(hoursx+1,1) = frequency(hoursx+1,1) + sum(p(timex,:,:)<1 & p(timex,:,:)>0,'all');
        frequency(hoursx+1,2) = frequency(hoursx+1,2) + sum(p(timex,:,:)<3 & p(timex,:,:)>=1, 'all');
        frequency(hoursx+1,3) = frequency(hoursx+1,3) + sum(p(timex,:,:)>=3, 'all');
    end
    
    
    for hoursx = hours
        % find hours wanted
        timex = hour(wtm)==hoursx;
        
        % create sum of precipitation for the year at that hour
        precipSum(hoursx+1, yrx-2000) = sum(p(timex,:,:),'all');
    end  
    
    clear p
    clear wtm
end


%% plot
disp('Plotting...');

hourly = sum(precipSum,2);
total = sum(precipSum,'all');
prcntSum = hourly ./ total;

fig1 = figure();
plot(0:23, (precipSum./sum(precipSum))');
xticks([0:2:23]);
xlim([0,23]);
xticklabels(["00","02","04","06","08","10","12","14","16","18","20","22"]);
xlabel("Hour (UTC)");
ylabel("Percentage of Total Precipitation (%)");
ylim([0,0.1]);
yticks([0,0.025,0.05,0.075,0.1]);
yticklabels(["0","2.5","5","7.5","10"]);
grid on
title("Total Yearly Precipitation at Each Hour");

fig_name = save_dir + "plot_prcntOfPrecip_" + string(yr_initial) + "_" + string(yr_final) + ".png";
export_fig(fig1, fig_name{:}, '-transparent', '-r300')


%% plot 

fig2 = figure();

ftotal = sum(frequency,'all');
prcnt = frequency ./ ftotal;

yyaxis left
p1 = plot(0:23,prcnt(:,1));
p1.LineWidth = 2;
xticks([0:2:23]);
xticklabels(["00","02","04","06","08","10","12","14","16","18","20","22"]);
xlim([0,23]);
ylabel("Frequency (%)",'Color', 'k');
yticks([0:0.01:0.08]);
yticklabels(["0","1","2","3","4","5","6","7","8"]);
xlabel("Hour (UTC)");
ylim([0,0.05]);
hold on
yyaxis right
p2 = plot(0:23,prcnt(:,2:3));
p2(1).LineWidth = 2;
p2(2).LineWidth = 2;
ylim([0 0.01]);
yticks([0,0.0025,0.005,0.0075,0.01]);
yticklabels(["0","0.25","0.5","0.75","1"]);

grid on
title("Frequency of Each Precipitation Event Type");
legend(["R < 1 mm/h","1 \leq R < 3 mm/h","R \geq 3 mm/h"]);

fig_name = save_dir + "plot_freqencyOfPrecip_" + string(yr_initial) + "_" + string(yr_final) + ".png";
export_fig(fig2, fig_name{:}, '-transparent', '-r300')



