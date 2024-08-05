%% Savanna Wolvin
% Created: Apr. 26th, 2023
% Edited: 

% SUMMARY

% INPUT

% OUTPUT


%% paths

clear;

% mapping functions
addpath('/uufs/chpc.utah.edu/common/home/strong-group4/savanna/m_map');

% export figure path
addpath('/uufs/chpc.utah.edu/common/home/strong-group4/savanna/export');

% countries shapefiles
addpath('/uufs/chpc.utah.edu/common/home/strong-group4/savanna/countries_shapefile');

% export figure path
addpath('/uufs/chpc.utah.edu/common/home/u1324060/nclColormap/');



%% variable presets

global var_name nice_name save_dir wrf_nest fnum prThrshld


corr_dir = "/uufs/chpc.utah.edu/common/home/strong-group4/savanna/himat_opg_poly/";
corr_file = "three_dim_lin_regress_";

wrf_nest = "wrfout_d02";

% radius value
radius = "25";

% facet data 
facet_dir   = "/uufs/chpc.utah.edu/common/home/strong-group4/savanna/wrf_facet_data/";

var_name = "P1";
nice_name = "p1";
var_num = [1,2];

% precip data
precip_dir  = "/uufs/chpc.utah.edu/common/home/strong-group4/savanna/precip_wrf/";

poly_fit = "polyfit2";

save_dir = "/uufs/chpc.utah.edu/common/home/strong-group4/savanna/himat_opg_poly/";


% Himalayas: 939, 955, 420, 795
% Karakoram: 388, 924, 526, 797
% Tibetan P: 975, 964, 423, 82
% Hindu-Kush: 12, 166, 501, 24
fnum = 939; % 955, 924, 964, 166

years = 2001:2015;

% precip lower threshold
prThrshld = 1;


%% map presets
disp('Load map data...');

global lati_d03 loni_d03 

% wrf file
fdir    = "/uufs/chpc.utah.edu/common/home/strong-group6/HIMAT/WRF_output/2001/output/";

% wrf file
fname   = "wrfout_d03_2001-01-01_00:00:00"; 
wrf_file    = fdir + fname;

lati_d03 = ncread(wrf_file,'XLAT',[1 1 1],[inf inf 1]); lati_d03 = double(lati_d03); 
loni_d03 = ncread(wrf_file,'XLONG',[1 1 1],[inf inf 1]); loni_d03 = double(loni_d03); 
elev = ncread(wrf_file,'HGT',[1 1 1],[inf inf 1]); elev = elev';


%% Load data

disp("Loading Facets, Polynomial Coefficient, & Atmospheric Predictor...");

load(facet_dir + 'wrf_facets_' + radius + '.mat');
global facet_contour days
use = (facets_labeled==fnum);
facet_contour = double(use);

poly_var = [];
days = [];
precip = [];
pr_daily = [];

for yr = years
    % load in precip
    precip_fname = "daily_precip_" + string(yr) + ".mat";
    precip_file = precip_dir + precip_fname;
    load(precip_file);
    precip = cat(1, precip, nanmean(pr(:,use==1), 2));
    pr_daily = cat(1, pr_daily, pr(:, use==1));


    load(corr_dir + radius + "km/daily_opg_" + string(yr) + "_" + radius + "km_fitnlm_nov19_bottomElevationZero.mat"); 
    varx = NaN(size(t,1),3);
    numDays = size(t,1);
    for dayx = 1:numDays
        varx(dayx,1) = mode(opg_polyfit2(1,dayx,use)); % p1
        varx(dayx,2) = mode(opg_polyfit2(2,dayx,use)); % p2
        varx(dayx,3) = mode(opg_polyfit2(3,dayx,use)); % p3
    end

    poly_var = cat(1,poly_var,varx);
    days = cat(1, days, t);

end


%% seasonal vectors
disp("Create Seasonal Vectors...");

djf_days = zeros(size(days,1),1);
mam_days = zeros(size(days,1),1);
jja_days = zeros(size(days,1),1);
son_days = zeros(size(days,1),1);

for day = 1:size(days,1)
    d_vec = datevec(days(day,1));
    if d_vec(1,2) == 1 || d_vec(1,2) == 2 || d_vec(1,2) == 12
        djf_days(day,1) = 1;
    elseif d_vec(1,2) == 3 || d_vec(1,2) == 4 || d_vec(1,2) == 5
        mam_days(day,1) = 1;
    elseif d_vec(1,2) == 6 || d_vec(1,2) == 7 || d_vec(1,2) == 8
        jja_days(day,1) = 1;
    elseif d_vec(1,2) == 9 || d_vec(1,2) == 10 || d_vec(1,2) == 11
        son_days(day,1) = 1;
    end
end

djf_days = logical(djf_days);
mam_days = logical(mam_days);
jja_days = logical(jja_days);
son_days = logical(son_days);
mam_son_days = logical(mam_days + son_days);


%% pull percentile values

p1_q33 = prctile(poly_var(:,1), (1/3)*100);
p1_q66 = prctile(poly_var(:,1), (2/3)*100);

p2_q33 = prctile(poly_var(:,2), (1/3)*100);
p2_q66 = prctile(poly_var(:,2), (2/3)*100);

p1 = poly_var(:,1);
p2 = poly_var(:,2);


%% Individual Case

%%% DJF
opg_djf = NaN(4,3);

% seasonal means
select_logical = ~isnan(poly_var(:,1)) & (djf_days) &...
    (precip(:) > prThrshld);
opg_djf(1,1) = mean(poly_var(select_logical,1));
opg_djf(1,2) = mean(poly_var(select_logical,2));
opg_djf(1,3) = mean(poly_var(select_logical,3));

% CASE 1
select_logical = ~isnan(poly_var(:,1)) & (djf_days) &...
    (precip(:) > prThrshld) & (poly_var(:,1) < p1_q33) & (poly_var(:,2) > p2_q66);
opg_djf(2,1) = mean(poly_var(select_logical,1));
opg_djf(2,2) = mean(poly_var(select_logical,2));
opg_djf(2,3) = mean(poly_var(select_logical,3));

% CASE 5
select_logical = ~isnan(poly_var(:,1)) & (djf_days) &...
    (precip(:) > prThrshld) & (poly_var(:,1) > p1_q33) & (poly_var(:,1) < p1_q66) &...
    (poly_var(:,2) > p2_q33) & (poly_var(:,2) < p2_q66);
opg_djf(3,1) = mean(poly_var(select_logical,1));
opg_djf(3,2) = mean(poly_var(select_logical,2));
opg_djf(3,3) = mean(poly_var(select_logical,3));

% CASE 9
select_logical = ~isnan(poly_var(:,1)) & (djf_days) &...
    (precip(:) > prThrshld) & (poly_var(:,1) > p1_q66) & (poly_var(:,2) < p2_q33);
opg_djf(4,1) = mean(poly_var(select_logical,1));
opg_djf(4,2) = mean(poly_var(select_logical,2));
opg_djf(4,3) = mean(poly_var(select_logical,3));

%%% JJA
opg_jja = NaN(4,3);

% seasonal means
select_logical = ~isnan(poly_var(:,1)) & (jja_days) &...
    (precip(:) > prThrshld);
opg_jja(1,1) = mean(poly_var(select_logical,1));
opg_jja(1,2) = mean(poly_var(select_logical,2));
opg_jja(1,3) = mean(poly_var(select_logical,3));

% CASE 1
select_logical = ~isnan(poly_var(:,1)) & (jja_days) &...
    (precip(:) > prThrshld) & (poly_var(:,1) < p1_q33) & (poly_var(:,2) > p2_q66);
opg_jja(2,1) = mean(poly_var(select_logical,1));
opg_jja(2,2) = mean(poly_var(select_logical,2));
opg_jja(2,3) = mean(poly_var(select_logical,3));

% CASE 5
select_logical = ~isnan(poly_var(:,1)) & (jja_days) &...
    (precip(:) > prThrshld) & (poly_var(:,1) > p1_q33) & (poly_var(:,1) < p1_q66) &...
    (poly_var(:,2) > p2_q33) & (poly_var(:,2) < p2_q66);
opg_jja(3,1) = mean(poly_var(select_logical,1));
opg_jja(3,2) = mean(poly_var(select_logical,2));
opg_jja(3,3) = mean(poly_var(select_logical,3));

% CASE 9
select_logical = ~isnan(poly_var(:,1)) & (jja_days) &...
    (precip(:) > prThrshld) & (poly_var(:,1) > p1_q66) & (poly_var(:,2) < p2_q33);
opg_jja(4,1) = mean(poly_var(select_logical,1));
opg_jja(4,2) = mean(poly_var(select_logical,2));
opg_jja(4,3) = mean(poly_var(select_logical,3));


% cmap for scatter
mamson_c = 'k';
djf_c    = [0.76, 0.227, 0];
jja_c    = [0.0941, 0.388, 0];


%% Combo plot, scatter, boxplot, average OPG

use = (facets_labeled == fnum);
x = elev(use);
x1 = [min(x):max(x)];
x1v2 = (x1-min(x))./1000;

y_bw = zeros(size(p1));
y_bw(mam_son_days) = 1;
y_bw(jja_days) = 0;
y_bw(djf_days) = -1;


% font size
ftsz = 8;

% tercile linewidth
lw = 1.5;

% scatter min/max
p1_lim = [-8,4];
p2_lim = [-40, 40];

% fig = figure('Visible', 'off');
fig = figure('Position', [10, 10, 400, 500]);

tt = tiledlayout(8,4);
tt.TileSpacing = 'none';


% Box plot p2
nexttile([3,1]);

patch([-1.75,1.75,1.75,-1.75],[p2_q33, p2_q33, p2_lim(1), p2_lim(1)],...
    'r', 'facealpha', 0, 'EdgeColor', 'b', 'LineWidth', lw);
patch([-1.75,1.75,1.75,-1.75],[p2_q33+.01, p2_q33+.01, p2_q66-.01, p2_q66-.01],...
    [0.66 0 1], 'facealpha', 0, 'EdgeColor', [0.66 0 1], 'LineWidth', lw);
patch([-1.75,1.75,1.75,-1.75],[p2_q66, p2_q66, p2_lim(2), p2_lim(2)],...
    'b', 'facealpha', 0, 'EdgeColor', 'r', 'LineWidth', lw);

hold on;

boxchart(y_bw, p2, 'BoxFaceColor','k' ,'Orientation','vertical', ...
    'JitterOutliers', 'on', 'MarkerStyle', '.', 'BoxWidth', 0.95,...
    'MarkerColor', [0.2,0.2,0.2], 'MarkerSize',4);

xlim([-1.75,1.75]);
ylim(p2_lim)
xticks([-1,0,1]);
xticklabels({'DJF','JJA','MAMSON'});

ylabel("OPG at the Lowest Elevation (P2)");

ax = gca; ax.FontSize = ftsz;
grid on;

% scatter plot
ax_sca = nexttile([3,3]);
scatter(p1, p2, 5, [0.2,0.2,0.2], 'filled', 'DisplayName', 'Obs');
hold on;

plot([min(p1), p1_q33, p1_q33], [p2_q66, p2_q66, max(p2)],...
    'r', 'LineWidth', lw, 'DisplayName', 'P1(Lower Tercile) & P2(Upper Tercile)');
plot([p1_q33, p1_q66, p1_q66, p1_q33, p1_q33], [p2_q33, p2_q33, p2_q66, p2_q66, p2_q33],...
    'color', [0.66, 0, 1], 'LineWidth', lw, 'DisplayName', 'P1(Middle Tercile) & P2(Middle Tercile)');
plot([max(p1), p1_q66, p1_q66], [p2_q33, p2_q33, min(p2)],...
    'b', 'LineWidth', lw, 'DisplayName', 'P1(Upper Tercile) & P2(Lower Tercile)');

% limits 
xlim(p1_lim);
ylim(p2_lim);

ax_sca.XTickLabel = [];
ax_sca.YTickLabel = [];

title("a) Relationship Between Facet-939 Coefficients");

grid on;

% skip tile 
ax_blank = nexttile([1,1]);
ax_blank.Visible = 'off';

% box plot p1
nexttile([1,3]);

patch([p1_q33, p1_q33, p1_lim(1), p1_lim(1)],[-1.75,1.75,1.75,-1.75], ...
    'r', 'facealpha', 0, 'EdgeColor', 'r', 'LineWidth', lw);
patch([p1_q33+.01, p1_q33+.01, p1_q66-.01, p1_q66-.01],[-1.75,1.75,1.75,-1.75], ...
    [0.66 0 1], 'facealpha', 0, 'EdgeColor', [0.66 0 1], 'LineWidth', lw);
patch([p1_q66, p1_q66, p1_lim(2), p1_lim(2)],[-1.75,1.75,1.75,-1.75], ...
    'b', 'facealpha', 0, 'EdgeColor', 'b', 'LineWidth', lw);

hold on;
boxchart(y_bw, p1, 'BoxFaceColor','k','Orientation','horizontal', ...
    'JitterOutliers', 'on', 'MarkerStyle', '.', 'BoxWidth', 0.95,...
    'MarkerColor', [0.2,0.2,0.2], 'MarkerSize',4);

xlim(p1_lim);
ylim([-1.75,1.75]);
yticks([-1,0,1]);
yticklabels({'DJF','JJA',''});

xlabel("OPG Concavity (P1)");

ax = gca; ax.FontSize = ftsz;
grid on;

% skip tile 
ax_blank = nexttile([1,4]);
ax_blank.Visible = 'off';

% facet 939 DJF
window = 200;

nexttile([2,2]);
%%% CASE 1
pp1 = plot(x1/1000, opg_djf(2,1)*(x1v2.^2) + opg_djf(2,2)*(x1v2) + opg_djf(2,3),'color', 'r');
pp1.Color(4) = 0.9;
pp1.LineWidth = 2;

hold on
select_logical = ~isnan(poly_var(:,1)) & (djf_days) &...
    (precip(:) > prThrshld) & (poly_var(:,1) < p1_q33) & (poly_var(:,2) > p2_q66);
pr_case1 = mean(pr_daily(select_logical,:), 1 , 'omitmissing');
[elev_sort, idx] = sort(elev(use)/1000);

mpre = movmean(pr_case1(idx), window);
plot(elev_sort, mpre, 'r--', 'LineWidth', 1);


%%% CASE 5
pp2 = plot(x1/1000, opg_djf(3,1)*(x1v2.^2) + opg_djf(3,2)*(x1v2) + opg_djf(3,3),'color', [0.66 0 1]);
pp2.Color(4) = 0.9;
pp2.LineWidth = 2;

select_logical = ~isnan(poly_var(:,1)) & (djf_days) &...
    (precip(:) > prThrshld) & (poly_var(:,1) > p1_q33) & (poly_var(:,1) < p1_q66) &...
    (poly_var(:,2) > p2_q33) & (poly_var(:,2) < p2_q66);
pr_case5 = mean(pr_daily(select_logical,:), 1 , 'omitmissing');

mpre = movmean(pr_case5(idx), window);
plot(elev_sort, mpre, 'Color', [0.66 0 1], 'LineStyle', '--', 'LineWidth', 1);

%%% CASE 9
pp3 = plot(x1/1000, opg_djf(4,1)*(x1v2.^2) + opg_djf(4,2)*(x1v2) + opg_djf(4,3),'color', 'b');
pp3.Color(4) = 0.9;
pp3.LineWidth = 2;

select_logical = ~isnan(poly_var(:,1)) & (djf_days) &...
    (precip(:) > prThrshld) & (poly_var(:,1) > p1_q66) & (poly_var(:,2) < p2_q33);
pr_case9 = mean(pr_daily(select_logical,:), 1 , 'omitmissing');

mpre = movmean(pr_case9(idx), window);
plot(elev_sort, mpre, 'b--', 'LineWidth', 1);

title({"c) Facet-" + num2str(fnum) + ": DJF"});

xlabel("Elevation (km)");
% ylim([0,40]);
xlim([min(elev(use)), max(elev(use))]/1000);
ylabel("Average Precip. Total (mm)");


% custom legend addition
plot(nan, 'k--', 'LineWidth', 1);

lgd3 = legend(["Strong Sublinear Increase", "",...
    "Sublinear Increase", "", ...
    "Superlinear Increase or Sublinear Decrease", "", "WRF Simulation Moving Average"], ...
    'location', 'southoutside', ...
    NumColumns=1, FontSize=ftsz);
legend('boxoff')

ax = gca; ax.FontSize = ftsz;

grid on;
hold off;

% facet 939 JJA
nexttile([2,2]);
%%% CASE 1
pp1 = plot(x1/1000, opg_jja(2,1)*(x1v2.^2) + opg_jja(2,2)*(x1v2) + opg_jja(2,3),'color', 'r');
pp1.Color(4) = 0.9;
pp1.LineWidth = 2;

hold on

select_logical = ~isnan(poly_var(:,1)) & (jja_days) &...
    (precip(:) > prThrshld) & (poly_var(:,1) < p1_q33) & (poly_var(:,2) > p2_q66);
pr_case1 = mean(pr_daily(select_logical,:), 1 , 'omitmissing');

mpre = movmean(pr_case1(idx), window);
plot(elev_sort, mpre, 'r--', 'LineWidth', 1);

%%% CASE 5
pp2 = plot(x1/1000, opg_jja(3,1)*(x1v2.^2) + opg_jja(3,2)*(x1v2) + opg_jja(3,3),'color', [0.66 0 1]);
pp2.Color(4) = 0.9;
pp2.LineWidth = 2;

select_logical = ~isnan(poly_var(:,1)) & (jja_days) &...
    (precip(:) > prThrshld) & (poly_var(:,1) > p1_q33) & (poly_var(:,1) < p1_q66) &...
    (poly_var(:,2) > p2_q33) & (poly_var(:,2) < p2_q66);
pr_case5 = mean(pr_daily(select_logical,:), 1 , 'omitmissing');

mpre = movmean(pr_case5(idx), window);
plot(elev_sort, mpre, 'Color', [0.66 0 1], 'LineStyle', '--', 'LineWidth', 1);


%%% CASE 3
pp3 = plot(x1/1000, opg_jja(4,1)*(x1v2.^2) + opg_jja(4,2)*(x1v2) + opg_jja(4,3),'color', 'b');
pp3.Color(4) = 0.9;
pp3.LineWidth = 2;

select_logical = ~isnan(poly_var(:,1)) & (jja_days) &...
    (precip(:) > prThrshld) & (poly_var(:,1) > p1_q66) & (poly_var(:,2) < p2_q33);
pr_case9 = mean(pr_daily(select_logical,:), 1 , 'omitmissing');

mpre = movmean(pr_case9(idx), window);
plot(elev_sort, mpre, 'b--', 'LineWidth', 1);

title({"d) Facet-" + num2str(fnum) + ": JJA"});

xlabel("Elevation (km)");
% ylim([0,40]);
yticklabels([]);
xlim([min(elev(use)), max(elev(use))]/1000);


ax = gca; ax.FontSize = ftsz;

grid on;
hold off;

% skip tile 
ax_blank = nexttile([1,4]);
ax_blank.Visible = 'off';





