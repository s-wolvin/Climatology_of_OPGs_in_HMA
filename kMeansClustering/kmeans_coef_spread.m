%% Savanna Wolvin
% Created: Jul. 8th, 2024
% Edited: Oct. 6th, 2022

% SUMMARY
% Plot spread of cluster coefficients by clustering the coefficients 1000
% times with seed values ranging from 1 to 1000. To match clusters of
% differing seeds, cluster centroids are compared


%% Add Paths

addpath('/uufs/chpc.utah.edu/common/home/strong-group4/savanna/m_map');
addpath('/uufs/chpc.utah.edu/common/home/strong-group4/savanna/countries_shapefile');
addpath('/uufs/chpc.utah.edu/common/home/u1324060/exp_fig');
addpath('/uufs/chpc.utah.edu/common/home/u1324060/nclColormap/');


%% User-defined values

max_cluster     = 10;

years           = 2001:2015;

radius          = "25";
poly_fit        = "fitnlm";
grad            = "opg";
grad_dir        = "/uufs/chpc.utah.edu/common/home/strong-group4/savanna/himat_opg_poly/" + radius + "km/";
poly_coef       = 1:3;
value           = "z-score_timeseries";
poly_name       = "z-score_timeseries";
clust_dir   = "/uufs/chpc.utah.edu/common/home/strong-group4/savanna/cluster/clustStruct/";
clust_dataset = "k-means_opg_coef_" + poly_name + "_value_" + value + "_" + radius + "km_fitnlm_nov19_bottomElevationZero_mean.mat";

% raw, z-score_timeseries, z-score_facets, z-score_all, normalizeX-Axis,
% median_timeseries
value           = "z-score_timeseries";



% facet data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
facet_dir   = "/uufs/chpc.utah.edu/common/home/strong-group4/savanna/wrf_facet_data/";


%% Leftover vars

months = [1:12];

%% BIG LOOOOP

% file folders %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig_folder      = ['/uufs/chpc.utah.edu/common/home/strong-group4/savanna/cluster/clustFigs/'];          % Path to folder to save figures
data_folder      = ['/uufs/chpc.utah.edu/common/home/strong-group4/savanna/cluster/clustStruct/']; 

% load in maps presets
fdir    = "/uufs/chpc.utah.edu/common/home/strong-group6/HIMAT/WRF_output/2000/output/";
fname3   = "wrfout_d03_2000-01-01_00:00:00"; 
wrf_file3    = fdir + fname3;

lat = ncread(wrf_file3,'XLAT',[1 1 1],[inf inf 1]); lat = double(lat'); 
lon = ncread(wrf_file3,'XLONG',[1 1 1],[inf inf 1]); lon = double(lon'); 

map_latlims     = [min(lat(:)) max(lat(:))]; 
map_lonlims     = [min(lon(:)) max(lon(:))];


% reshape to make the locations be listed in rows and not a 3D matrix of
% Facets
facet_fname = "wrf_facets_" + radius;
load(facet_dir + facet_fname);

size_fi = NaN(1,max(facets_labeled(:)));
for fi = 1:max(facets_labeled)
    size_fi(fi) = sum(facets_labeled(:)==fi);
end
size_fi = size_fi';

opg_opg_monthly = [];

for coef_num = poly_coef
    opg_polyx = [];
    d_vec = [];
    
    disp('Load Coefficient Number ' + string(coef_num))
    for yr = years
        poly_values = "daily_" + grad + "_" + yr + "_" + radius + "km_" + poly_fit + "_nov19_bottomElevationZero.mat";
        load(grad_dir + poly_values);
        opg_polyx = cat(2, opg_polyx, opg_polyfit2(coef_num, :, :, :));

        d_vec = cat(1, d_vec, datevec(t));
    end
    
    disp("Reshape Matrix into 2D");
    opg_polyx = squeeze(opg_polyx);
    polyVar_2D = nan(max(facets_labeled(:)) , size(opg_polyx, 1));
    for fi = 1:max(facets_labeled(:)) 
        use = (facets_labeled == fi);
        rowA = opg_polyx(:, use);
        polyVar_2D(fi, :) = rowA(:,1);
    end

    % Create Average Field
    disp('Monthly Mean...');

    polyVar_2D_monthlyx = zeros(size(polyVar_2D, 1), length(months));

    if value == "z-score_timeseries"
        polySTD = std(polyVar_2D, 0, 2, 'omitnan');
        polyMean = mean(polyVar_2D, 2, 'omitnan');
        polyVar_2D = (polyVar_2D - polyMean) ./ polySTD;
    elseif value == "z-score_facets"
        polySTD = std(polyVar_2D, 0, 1, 'omitnan');
        polyMean = mean(polyVar_2D, 1, 'omitnan');
        polyVar_2D = (polyVar_2D - polyMean) ./ polySTD;
    elseif value == "z-score_all"
        polySTD = std(polyVar_2D(:), 'omitnan');
        polyMean = mean(polyVar_2D(:), 'omitnan');
        polyVar_2D = (polyVar_2D - polyMean) ./ polySTD;
    elseif value == "normalizeX-Axis"
        fdir    = "/uufs/chpc.utah.edu/common/home/strong-group6/HIMAT/WRF_output/2000/output/";
        fname   = "wrfout_d03_2000-01-01_00:00:00"; 
        elev = ncread(fdir + fname, 'HGT', [1 1 1], [inf inf 1])';
        elev = elev ./ 1000;

        if coef_num == 1
            for fi = 1:size(polyVar_2D, 1)
                use = (facets_labeled == fi);
                dElev = max(elev(use)) - min(elev(use));
                polyVar_2D(fi, :) = polyVar_2D(fi, :) .* dElev .* dElev;
            end
        elseif coef_num == 2
            for fi = 1:size(polyVar_2D, 1)
                use = (facets_labeled == fi);
                dElev = max(elev(use)) - min(elev(use));
                polyVar_2D(fi, :) = polyVar_2D(fi, :) .* dElev;
            end
        end

    end

    for fi = 1:max(facets_labeled(:))
        count = 1;
        for mx = months
            polyVar_2D_monthlyx(fi, count) =  mean(polyVar_2D(fi, d_vec(:,2)==mx), 'omitnan');
            % polyVar_2D_monthlyx(fi, count) =  median(polyVar_2D(fi, d_vec(:,2)==mx), 'omitnan');
            count = count + 1;
        end
    end

    opg_opg_monthly = cat(2, opg_opg_monthly, polyVar_2D_monthlyx);
end


%%

C = [rgb('cornflowerblue');...
    rgb('navajowhite');...
    rgb('darkorchid');...
    rgb('deeppink');... 
    rgb('lightseagreen');...
    rgb('teal');...
    rgb('firebrick');];


%% create spread of clustering

% load in clusters
load(clust_dir + clust_dataset);
% clusters = cluster_data.clusters;
centroids_orig = cluster_data.centroids(5,1);
centroids_orig = centroids_orig{1,1};
clust_layer = cluster_data.clusters(:,5);


cnum = 5;
seed_num = 1000;
clust_spread = NaN(size(opg_opg_monthly,1), seed_num);

% rng(1);
for seed = 1:seed_num
    % disp(seed);
    orig = [1,2,3,4,5];
    new = [NaN, NaN, NaN, NaN, NaN];
    clust_idx = [true,true,true,true,true];

    rng(seed);
    [idx, centroid, sumD] = kmeans(opg_opg_monthly, cnum);


    %%%%%%%%%%%%%%%%%%%%%% Matching clusters based on facets
    % for clust_orig = 1:5
    %     log_orig = [(clust_layer == 1), (clust_layer == 2), (clust_layer == 3), (clust_layer == 4), (clust_layer == 5)];
    %     log_new = [(idx == 1), (idx == 2), (idx == 3), (idx == 4), (idx == 5)];
    % 
    %     % formulate difference
    %     [~, dist] = dsearchn(log_orig(:,clust_orig)', log_new(:,clust_idx)');
    % 
    %     [m, i] = min(dist);
    %     where = orig(clust_idx);
    %     new(:,clust_orig) = where(i);
    %     clust_idx(:,where(i)) = false;
    % 
    % end
    % idx_new = NaN(size(idx));
    % idx_new(idx==new(1)) = 1;
    % idx_new(idx==new(2)) = 2;
    % idx_new(idx==new(3)) = 3;
    % idx_new(idx==new(4)) = 4;
    % idx_new(idx==new(5)) = 5;
    %%%%%%%%%%%%%%%%%%
    

    %%%%%%%%%%%%%%%%%%%%% Matching clusters based on the centroid
    for clust_orig = 1:5
        [~, dist] = dsearchn(centroids_orig(clust_orig,:), centroid(clust_idx,:));
        [~, i] = min(dist);
        where = orig(clust_idx);
        new(:,clust_orig) = where(i);
        clust_idx(:,where(i)) = false;
    end
    idx_new = NaN(size(idx));
    idx_new(idx==new(1)) = 1;
    idx_new(idx==new(2)) = 2;
    idx_new(idx==new(3)) = 3;
    idx_new(idx==new(4)) = 4;
    idx_new(idx==new(5)) = 5;

    clust_spread(:, seed) = idx_new;
end







%% load coef values


% reshape to make the locations be listed in rows and not a 3D matrix of
% Facets
% facet_fname = "wrf_facets_" + radius;
% load(facet_dir + facet_fname);


opg_opg_monthly_plot = [];

for coef_num = poly_coef
    opg_polyx = [];
    d_vec = [];

    disp('Load Coefficient Number ' + string(coef_num))
    for yr = years
        poly_values = "daily_" + grad + "_" + yr + "_" + radius + "km_" + poly_fit + "_nov19_bottomElevationZero.mat";
        load(grad_dir + poly_values);
        opg_polyx = cat(2, opg_polyx, opg_polyfit2(coef_num, :, :, :));

        d_vec = cat(1, d_vec, datevec(t));
    end

    disp("Reshape Matrix into 2D");
    opg_polyx = squeeze(opg_polyx);
    polyVar_2D = nan(max(facets_labeled(:)) , size(opg_polyx, 1));
    for fi = 1:max(facets_labeled(:)) 
        use = (facets_labeled == fi);
        rowA = opg_polyx(:, use);
        polyVar_2D(fi, :) = rowA(:,1);
    end

    % Create Average Field
    disp('Monthly Mean...');

    polyVar_2D_monthlyx = zeros(size(polyVar_2D, 1), length(months));

    for fi = 1:max(facets_labeled(:))
        count = 1;
        for mx = months
            polyVar_2D_monthlyx(fi, count) =  mean(polyVar_2D(fi, d_vec(:,2)==mx), 'omitnan');
            count = count + 1;
        end
    end

    opg_opg_monthly_plot = cat(2, opg_opg_monthly_plot, polyVar_2D_monthlyx);
end


%% tiled plot of the coefficients MEDIAN

fig1 = figure('Position',[10 10 550 650]);

tt = tiledlayout(3,1);
tt.TileSpacing = 'compact';
tt.Padding = 'compact';

% Tile 1
nexttile

opg_median_range = NaN(seed_num, 12);

for cx = 1:cnum
    for seed = 1:seed_num
        clust_seed = clust_spread(:, seed);
        opg_median_range(seed, :) = median(opg_opg_monthly_plot(clust_seed==cx, 1:12), 'omitmissing');
    end

    plot(1:12, median(opg_median_range), color=faceColor(cx), LineWidth=3);

    hold on;

    plot(1:12, prctile(opg_median_range, 5), "--", color=faceColor(cx), LineWidth=1);
    plot(1:12, prctile(opg_median_range, 95), "--", color=faceColor(cx), LineWidth=1);
end

xlim([1,12]);
xticks([1 2 3 4 5 6 7 8 9 10 11 12]);
set(gca, 'XTickLabel', {'J' 'F' 'M' 'A' 'M' 'J' 'J' 'A' 'S' 'O' 'N' 'D'}, ...
    'fontsize', 11);
ylim([-2.5,2.5]);
title("a) OPG Concavity (P1)");
ylabel("mm/m^2");
grid on;

% Tile 2
nexttile

opg_median_range = NaN(seed_num, 12);

for cx = 1:cnum
    for seed = 1:seed_num
        clust_seed = clust_spread(:, seed);
        opg_median_range(seed, :) = median(opg_opg_monthly_plot(clust_seed==cx, 13:24), 'omitmissing');
    end

    plot(1:12, median(opg_median_range), color=faceColor(cx), LineWidth=3);

    hold on;

    plot(1:12, prctile(opg_median_range, 5), "--", color=faceColor(cx), LineWidth=1);
    plot(1:12, prctile(opg_median_range, 95), "--", color=faceColor(cx), LineWidth=1);
end

xlim([1,12]);
xticks([1 2 3 4 5 6 7 8 9 10 11 12]);
set(gca, 'XTickLabel', {'J' 'F' 'M' 'A' 'M' 'J' 'J' 'A' 'S' 'O' 'N' 'D'}, ...
    'fontsize', 11);
ylim([-0.75,7.25]);
title("b) OPG at the Facet's Lowest Elevation (P2)");
ylabel("mm/m");
grid on;



% Tile 3
nexttile

opg_median_range = NaN(seed_num, 12);

for cx = 1:cnum
    for seed = 1:seed_num
        clust_seed = clust_spread(:, seed);
        opg_median_range(seed, :) = median(opg_opg_monthly_plot(clust_seed==cx, 25:36), 'omitmissing');
    end

    plot(1:12, median(opg_median_range), color=faceColor(cx), LineWidth=3);

    hold on;

    plot(1:12, prctile(opg_median_range, 5), "--", color=faceColor(cx), LineWidth=1);
    plot(1:12, prctile(opg_median_range, 95), "--", color=faceColor(cx), LineWidth=1);
end

xlim([1,12]);
xticks([1 2 3 4 5 6 7 8 9 10 11 12]);
set(gca, 'XTickLabel', {'J' 'F' 'M' 'A' 'M' 'J' 'J' 'A' 'S' 'O' 'N' 'D'}, ...
    'fontsize', 11);
ylim([0,5.5]);
title("c) Precip. Total at the Facet's Lowest Elevation (P3)");
ylabel("mm");
grid on;

title(tt, "Median of Monthly OPG Coefficients");

% create a legend
qw{1} = plot(nan, color=faceColor(1), LineWidth=3);
qw{2} = plot(nan, color=faceColor(2), LineWidth=3);
qw{3} = plot(nan, color=faceColor(3), LineWidth=3);
qw{4} = plot(nan, color=faceColor(4), LineWidth=3); 
qw{5} = plot(nan, color=faceColor(5), LineWidth=3); 
qw{6} = plot(nan, 'k--'); 
legend([qw{:}], {'Cluster 1','Cluster 2','Cluster 3', 'Cluster 4', 'Cluster 5', '5th to 95th Percentile'}, ...
    'location', 'southoutside', 'NumColumns', 3)


% savefigname = ["/uufs/chpc.utah.edu/common/home/u1324060/himat_ms/pub_fig/" +...
%         grad + "_k-means_" + poly_name + "_for_" + cluster_number + "_median_coef.png"];
% export_fig(savefigname{:}, '-transparent', '-r500');



%% tiled plot of the coefficients of the median, 25th, and 75th percentile

fig1 = figure('Position',[10 10 600 650]);

tt = tiledlayout(3,5);
tt.TileSpacing = 'none';
tt.Padding = 'compact';

% Tile 1

opg_median_range = NaN(seed_num, 12);
tile = [1,2,3,4,5];
clust = {"Cluster 1","Cluster 2","Cluster 3","Cluster 4","Cluster 5"};
lbls = {'a)', 'b)', 'c)', 'd)', 'e)'};

for cx = 1:cnum
    for seed = 1:seed_num
        clust_seed = clust_spread(:, seed);
        opg_median_range(seed, :) = median(opg_opg_monthly_plot(clust_seed==cx, 1:12), 'omitmissing');
    end
    
    nexttile(tile(cx))

    patch([1:12, 12:-1:1], [prctile(opg_median_range, 5), flip(prctile(opg_median_range, 95))], ...
        faceColor(cx), 'EdgeColor', [0.5,0.5,0.5]);

    hold on;

    plot(1:12, median(opg_median_range), color='k', LineWidth=2);


    xlim([1,12]);
    xticks([1 2 3 4 5 6 7 8 9 10 11 12]);
    xticklabels([]);
    set(gca, 'fontsize', 9, 'YColor', 'k', 'XColor', 'k');
    ylim([-2.5,2.5]);
    % title(clust(cx), 'FontWeight','normal');
    text(6, 2.25, lbls(cx), 'fontsize', 11);
    grid on;

    if cx == 1
        ylabel("P1 (mm/m^2)");
    else
        yticklabels([]);
    end

end




% Tile 2

opg_median_range = NaN(seed_num, 12);
tile = [6,7,8,9,10];
lbls = {'f)', 'g)', 'h)', 'i)', 'j)'};

for cx = 1:cnum
    for seed = 1:seed_num
        clust_seed = clust_spread(:, seed);
        opg_median_range(seed, :) = median(opg_opg_monthly_plot(clust_seed==cx, 13:24), 'omitmissing');
    end

    nexttile(tile(cx))
    
    patch([1:12, 12:-1:1], [prctile(opg_median_range, 5), flip(prctile(opg_median_range, 95))], ...
        faceColor(cx), 'EdgeColor', [0.5,0.5,0.5]);

    hold on;

    plot(1:12, median(opg_median_range), color='k', LineWidth=2);


    xlim([1,12]);
    xticks([1 2 3 4 5 6 7 8 9 10 11 12]);
    xticklabels([]);
    set(gca, 'fontsize', 9, 'YColor', 'k', 'XColor', 'k');
    ylim([-0.75,7.25]);
    text(6, 6.75, lbls(cx), 'fontsize', 11);
    grid on;

    if cx == 1
        ylabel("P2 (mm/m)");
    else
        yticklabels([]);
    end
end



% Tile 3

opg_median_range = NaN(seed_num, 12);
tile = [11,12,13,14,15];
lbls = {'k)', 'l)', 'm)', 'n)', 'o)'};

for cx = 1:cnum
    for seed = 1:seed_num
        clust_seed = clust_spread(:, seed);
        opg_median_range(seed, :) = median(opg_opg_monthly_plot(clust_seed==cx, 25:36), 'omitmissing');
    end

    nexttile(tile(cx))
    
    patch([1:12, 12:-1:1], [prctile(opg_median_range, 5), flip(prctile(opg_median_range, 95))], ...
        faceColor(cx), 'EdgeColor', [0.5,0.5,0.5]);

    hold on;

    plot(1:12, median(opg_median_range), color='k', LineWidth=2);

    xlim([1,12]);
    xticks([1 2 3 4 5 6 7 8 9 10 11 12]);
    set(gca, 'fontsize', 9, 'YColor', 'k', 'XColor', 'k');
    xticklabels({'Jan', '', '', '', 'May', '', '', '', 'Sep', '', '', ''});
    xtickangle(0);
    ylim([0,5.5]);
    text(6, 5,lbls(cx), 'fontsize', 11);
    grid on;

    if cx == 1
        ylabel("P3 (mm)");
    elseif cx == 3
        yticklabels([]);

        % create a legend
        qw{1} = plot(nan, color=faceColor(1), LineWidth=5);
        qw{2} = plot(nan, color=faceColor(2), LineWidth=5);
        qw{3} = plot(nan, color=faceColor(3), LineWidth=5);
        qw{4} = plot(nan, color=faceColor(4), LineWidth=5); 
        qw{5} = plot(nan, color=faceColor(5), LineWidth=5); 
        qw{6} = plot(nan, 'k', LineWidth=3); 
        qw{7} = plot(nan, color=[0.5,0.5,0.5], LineWidth=3); 
        legend([qw{:}], {'Cluster 1','Cluster 2', 'Cluster 3', 'Cluster 4', ...
            'Cluster 5', 'Median', '5th-95th Percentile'}, ...
            'location', 'southoutside', 'NumColumns', 3, 'FontSize', 10)
    else
        yticklabels([]);
    end
end


title(tt, "Spread of the Median of Monthly OPG Coefficients");



%% Color list function

function cmap = faceColor(num)
    C = [rgb('cornflowerblue');...
        rgb('navajowhite');...
        rgb('darkorchid');...
        rgb('deeppink');... 
        rgb('lightseagreen');...
        rgb('teal');...
        rgb('firebrick');...
        rgb('darkslateblue');...
        rgb('gold');...
        rgb('mistyrose');...
        rgb('lightslategrey');...
        rgb('red');...
        rgb('lightsalmon');...
        rgb('darkred');...
        rgb('mediumvioletred');...
        rgb('coral');...
        rgb('darkorange');...
        rgb('maroon');...
        rgb('aqua');...
        rgb('darkcyan');...
        rgb('darkgreen');...
        rgb('olivedrab');...
        rgb('teal');...
        rgb('saddlebrown');...
        rgb('tan');...
        rgb('peachpuff');...
        rgb('rosybrown');...
        rgb('midnightblue')];

    cmap = C(num,:);
    return
end
