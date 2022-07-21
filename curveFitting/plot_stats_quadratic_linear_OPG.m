%% Savanna Wolvin
% Created: Nov 20th, 2021
% Edited: Jul 21st, 2022

% SUMMARY
% Plot the pvalues and compare the R^2 values

% INPUT
% years -       Years of precipitation data to loop through
% radii -       Radius or radii used for faceting the terrain from
%               wrf_facets_MAIN.m
% p_threshold - The average facet precipiation needed for it to be 
%               defined as a preciaption event
% stats_dir -   Directory of the OPG fits stats data
% facet_dir -   Directory of the facet data
% wrf_dir -     Directory of the WRF output data
% wrf_fname -   File name of one WRF output file
% save_dir -    Directory to save to

% OUTPUT
% The output are four tiled plots:
%   precip_pvalue.png - percentage of precipitation events with a
%                       significant p-value for the first coefficient of 
%                       the quadratic fit
%   precip_aicc.png -   percentage of precipitation events with a lower or
%                       equal AICc-quadratic to that of the AICc-linear
%   precip_rsquared.png - percentange of precipitation events with a higher
%                       R^2 value for quadratic fits versus linear fits
%   precip_rsquared_diff.png - average percent difference between the
%                       quadratic and linear R^2 values for all
%                       precipitation events

% All four of these tiled plots show percentages in repect to each facet
% plotted on a map (left), and as a histogram by facet counts and grid
% point counts (right)




%% Add Paths

% mapping functions
addpath('/uufs/chpc.utah.edu/common/home/strong-group4/savanna/m_map');

% export figure path
addpath('/uufs/chpc.utah.edu/common/home/u1324060/exp_fig');

% export figure path
addpath('/uufs/chpc.utah.edu/common/home/u1324060/nclColormap');

% countries shapefiles
addpath('/uufs/chpc.utah.edu/common/home/strong-group4/savanna/countries_shapefile');




%% preset values

% years of events to evaluate
years       = 2001:2015;

% sampling radius of the faceting algorithm
radii       = "25";

% mean facet precipitation threshold
p_threshold = 0.10;

% folder for OPG fits
stats_dir    = "/uufs/chpc.utah.edu/common/home/strong-group4/savanna/himat_opg_poly/";

% folder for facet data
facet_dir   = "/uufs/chpc.utah.edu/common/home/strong-group4/savanna/wrf_facet_data/";

% WRF folder and file for elevation, lat, lon data
wrf_dir    = "/uufs/chpc.utah.edu/common/home/strong-group6/HIMAT/WRF_output/2000/output/";
wrf_fname   = "wrfout_d03_2000-01-01_00:00:00"; 

% save folder
save_dir    = "/uufs/chpc.utah.edu/common/home/strong-group4/savanna/comp_lin_poly/";




%% map presets
disp('Load map and facet data...');

% wrf file
wrf_file    = wrf_dir + wrf_fname;

% load wrf
lati = ncread(wrf_file,'XLAT',[1 1 1],[inf inf 1]); lati = double(lati); 
loni = ncread(wrf_file,'XLONG',[1 1 1],[inf inf 1]); loni = double(loni); 
elev = ncread(wrf_file, 'HGT', [1 1 1], [inf inf 1]);
elev = elev';
latlim = [min(lati(:)) max(lati(:))];
lonlim = [min(loni(:)) max(loni(:))];

% load facets
facet_fname = "wrf_facets_" + radii;
load(facet_dir + facet_fname); facets = facets';

% Load shapefiles
M=m_shaperead('ne_50m_admin_0_countries');




%% Load Data
disp("Load P-Value, R^2, and AICc values...");

% empty variables
lin_aicc = [];
lin_rsqu = [];

pol_aicc = [];
pol_rsqu = [];
pol_pval = [];

% Loop through each year, loading in the P-Values, R^2 values, and AICc
% indicies
for yr = years
    line_file = stats_dir + radii + "km/daily_opg_" + num2str(yr) + "_" + radii + "km_fitlm_stats_nov19_bottomElevationZero.mat";
    load(line_file);
    lin_aicc = cat(1, lin_aicc, aicc_lin);
    lin_rsqu = cat(1, lin_rsqu, r_squared_lin);
    
    
    poly_file = stats_dir + radii + "km/daily_opg_" + num2str(yr) + "_" + radii + "km_fitnlm_stats_nov19_bottomElevationZero.mat";
    load(poly_file);
    pol_aicc = cat(1, pol_aicc, aicc);
    pol_rsqu = cat(1, pol_rsqu, r_squared);
    pol_pval = cat(1, pol_pval, squeeze(p_value(1,:,:)));

end

clear aicc_lin r_squared_lin aicc r_squared p_value




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EVALUATE P-VALUE: MAPS THE PERCENT OF PRECIP EVENTS THAT ARE SIGNIFICANT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate precentage and create matrix for plotting

% Create Binary Matrix
binary_pval = pol_pval;
is_significant  = binary_pval < p_threshold;
not_significant = binary_pval >=  p_threshold;

binary_pval(not_significant) = 0;
binary_pval(is_significant) = 1;

clear is_significant not_significant

% Calculate Percentage
percnt = sum(binary_pval, 1, 'omitnan');
percnt = percnt ./ sum(~isnan(binary_pval),1);

% Create plotting matrix
perct_matrix = NaN(size(facets_labeled));
for fi = 1:max(facets_labeled(:))
    use = facets_labeled == fi;
    perct_matrix(use) = percnt(fi);
end



 
%% Plot Percentages of precip events with significant P-Values

fig1 = figure('Position',[10 10 1300 500]); 
t = tiledlayout(1, 2);
t.TileSpacing = 'compact';

% Map Plot
nexttile;
m_proj('sinusoidal', 'lon', lonlim, 'lat', latlim);
m_coast('patch', [.7 .7 .7]);

m_pcolor(loni, lati, perct_matrix'); 
for k=1:length(M.ncst)
     m_line(M.ncst{k}(:,1),M.ncst{k}(:,2), 'color', 'k','LineWidth',1.5); 
end

caxis([0, 1]);
colormap(BlueDarkRed18('Reverse', 1, 'Total', 10));
colorbar('southoutside', 'XTick', [0, 0.2, .4, 0.6, 0.8 1],...
    'XTickLabel',["0%","20%","40%","60%","80%","100%"],'FontSize',12);

m_grid('xtick', [70 80], 'ytick', [24 30 36], 'yaxislocation', 'left',...
    'xaxislocation', 'bottom','backgroundcolor',[0.8 0.8 0.8],'fontsize', 12);

% Histogram
nexttile;
yyaxis left;
histogram(percnt, 0:0.1:1, 'FaceColor', [0.0666 0.0863 0.9686]);
ylabel("Count of Facets");
ylim([0,250]);
hold on;

yyaxis right;
histogram(perct_matrix(:), 0:0.1:1, 'FaceColor', [0.9804 0.4666 0.0980]);
ylabel("Count of Grid Points");
ylim([0, 70000]);

xlim([0,1]);
xlabel("Percentage of Events",'FontSize',12);
xticks(0:0.2:1);
xticklabels(["0%","20%","40%","60%","80%","100%"]);
set(gca,'FontSize',12);

grid on;

title(t, "Percentage of Leading Coefficients with a p-value < " + string(p_threshold),...
    'fontsize', 16)

fig_name = save_dir + "precip_pvalue.png";
export_fig(fig1, fig_name{:}, '-transparent', '-r300')




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EVALUATE AICC: MAP PERCENTAGE OF PRECIP EVENTS WHERE QUAD-AICC IS LOWEST
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate precentage and create matrix for plotting

% Bianry Matrix for When the Polynomial AICC is Lower than Linear
pol_low = double(pol_aicc <= lin_aicc);
pol_low(isnan(pol_aicc)) = NaN;

% Calculate Percentage
percnt = sum(pol_low, 1, 'omitnan');
percnt = percnt ./ sum(~isnan(pol_low),1);

% Create plotting matrix
perct_matrix = NaN(size(facets_labeled));
for fi = 1:max(facets_labeled(:))
    use = facets_labeled == fi;
    perct_matrix(use) = percnt(fi);
end




%% Plot percentage of precip events where quadratic AICc is lower or equal 
% to the linear AICc

fig2 = figure('Position',[10 10 1300 500]); 
t = tiledlayout(1, 2);
t.TileSpacing = 'compact';

% Map Plot
nexttile;
m_proj('sinusoidal', 'lon', lonlim, 'lat', latlim);
m_coast('patch', [.7 .7 .7]);

m_pcolor(loni, lati, perct_matrix'); 
for k=1:length(M.ncst)
     m_line(M.ncst{k}(:,1),M.ncst{k}(:,2), 'color', 'k','LineWidth',1.5); 
end

caxis([0, 1]);
colormap(BlueDarkRed18('Reverse', 1, 'Total', 10));
colorbar('southoutside', 'XTick', [0, 0.2, .4, 0.6, 0.8 1],...
    'XTickLabel',["0%","20%","40%","60%","80%","100%"],'FontSize',12);

m_grid('xtick', [70 80], 'ytick', [24 30 36], 'yaxislocation', 'left',...
    'xaxislocation', 'bottom','backgroundcolor',[0.8 0.8 0.8],'fontsize', 12);

% Histogram of Percentage
nexttile; 
yyaxis left;
histogram(percnt, 0:0.1:1, 'FaceColor', [0.0666 0.0863 0.9686]);
ylabel("Count of Facets");
ylim([0,250]);
hold on;

yyaxis right;
histogram(perct_matrix(:), 0:0.1:1, 'FaceColor', [0.9804 0.4666 0.0980]);
ylabel("Count of Grid Points");
ylim([0,70000]);

xlim([0,1]);
xlabel("Percentage");
xticks(0:0.2:1);
xticklabels(["0%","20%","40%","60%","80%","100%"]);
set(gca,'FontSize',12);

grid on;

title(t, "Percentage of Precipitation Events with a Quadratic AICc Value Lower/Equal to a Linear AICc Value",...
    'fontsize', 16)

fig_name = save_dir + "precip_aicc.png";
export_fig(fig2, fig_name{:}, '-transparent', '-r300')




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EVALUATE R-SQUARED: MAP PERCENT OF PRECIP EVENTS WHERE QUAD-R^2 IS HIGHER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate precentage and create matrix for plotting

% Create binary matrix
pol_high = double(pol_rsqu > lin_rsqu);
pol_high(isnan(pol_rsqu)) = NaN;

% Calculate Percentage
percnt = sum(pol_high, 1, 'omitnan');
percnt = percnt ./ sum(~isnan(pol_high),1);

% Create plotting matrix
perct_matrix = NaN(size(facets_labeled));
for fi = 1:max(facets_labeled(:))
    use = facets_labeled == fi;
    perct_matrix(use) = percnt(fi);
end




%% Plot percentage of events where the r-squared values for a quadratic fit 
% is larger than the linear fit

fig3 = figure('Position',[10 10 1300 500]); 
t = tiledlayout(1, 2);
t.TileSpacing = 'compact';

% Map Plot
nexttile;
m_proj('sinusoidal', 'lon', lonlim, 'lat', latlim);
m_coast('patch', [.7 .7 .7]);

m_pcolor(loni, lati, perct_matrix'); 
for k=1:length(M.ncst)
     m_line(M.ncst{k}(:,1),M.ncst{k}(:,2), 'color', 'k','LineWidth',1.5); 
end

caxis([0, 1]);
colormap(BlueDarkRed18('Reverse', 1, 'Total', 10));
colorbar('southoutside', 'XTick', [0, 0.2, .4, 0.6, 0.8 1],...
    'XTickLabel',["0%","20%","40%","60%","80%","100%"],'FontSize',12);

m_grid('xtick', [70 80], 'ytick', [24 30 36], 'yaxislocation', 'left',...
    'xaxislocation', 'bottom','backgroundcolor',[0.8 0.8 0.8],'fontsize', 12);

% Histogram of Percentage
nexttile;
yyaxis left;
histogram(percnt, 0:0.1:1, 'FaceColor', [0.0666 0.0863 0.9686]);
ylabel("Count of Facets");
ylim([0,700]);
hold on;

yyaxis right;
histogram(perct_matrix(:), 0:0.1:1, 'FaceColor', [0.9804 0.4666 0.0980]);
ylabel("Count of Grid Points");
ylim([0,140000]);

xlim([0,1]);
xlabel("Percentage");
xticks(0:0.2:1);
xticklabels(["0%","20%","40%","60%","80%","100%"]);
set(gca,'FontSize',12);

grid on;

title(t, "Percentage of Precipitation Events where the R^2 Value is Higher for a Quadratic Model",...
    'fontsize', 16)


fig_name = save_dir + "precip_rsquared.png";
export_fig(fig3, fig_name{:}, '-transparent', '-r300')





%% Calculate difference and create matrix for plotting

% Bianry Matrix for When the Polynomial AICC is Lower than Linear
pol_high = pol_rsqu - lin_rsqu;

% Calculate Mean R squared difference
mean_val = mean(pol_high, 1, 'omitnan');

% Create plotting matrix
mean_matrix = NaN(size(facets_labeled));
for fi = 1:max(facets_labeled(:))
    use = facets_labeled == fi;
    mean_matrix(use) = mean_val(fi);
end




%% Plot percentage

fig4 = figure('Position',[10 10 1300 500]); 
t = tiledlayout(1, 2);
t.TileSpacing = 'compact';

% Map Plot
ax1 = nexttile;
m_proj('sinusoidal', 'lon', lonlim, 'lat', latlim);
m_coast('patch', [.7 .7 .7]);

m_pcolor(loni, lati, mean_matrix'); 
for k=1:length(M.ncst)
     m_line(M.ncst{k}(:,1),M.ncst{k}(:,2), 'color', 'k','LineWidth',1.5); 
end

caxis([-0.25, 0.25]);
colormap(BlueDarkRed18('Reverse', 1, 'Total', 10));
colorbar('southoutside', 'XTick', [-0.2, -0.1, 0, 0.1, 0.2],...
    'XTickLabel',["-20%","-10%","0%","10%","20%"],'FontSize',12);

m_grid('xtick', [70 80], 'ytick', [24 30 36], 'yaxislocation', 'left',...
    'xaxislocation', 'bottom','backgroundcolor',[0.8 0.8 0.8],'fontsize', 12);

% Histogram of Differences
ax2 = nexttile;
yyaxis left;
histogram(mean_val, -0.5:0.05:0.5, 'FaceColor', [0.0666 0.0863 0.9686]);
ylabel("Count of Facets");
ylim([0,700]);
hold on;

yyaxis right;
histogram(mean_matrix(:), -0.5:0.05:0.5, 'FaceColor', [0.9804 0.4666 0.0980]);
ylabel("Count of Grid Points");
ylim([0,140000]);

xlim([-0.25,0.25]);
xlabel("Difference");
xticks(-0.5:0.1:0.5);
xticklabels(["-50%","-40%","-30%","-20%", "-10%", "0%","10%","20%","30%","40%","50%"]);
set(gca,'FontSize',12);

grid on;

title(t, "Average Difference Between Quadratic and Linear R^2 Values of Precipitation Events",...
    'fontsize', 16)

fig_name = save_dir + "precip_rsquared_diff.png";
export_fig(fig4, fig_name{:}, '-transparent', '-r300')



