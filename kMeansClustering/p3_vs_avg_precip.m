%% Savanna Wolvin

% SUMMARY
% Plot the median of each cluster's P3 coefficient along with the mean
% precipitation total for each cluster. Additionally print the correlation
% between these two values


%% Variable prests

value           = "z-score_timeseries";
radius          = "25";
facet_dir       = "/uufs/chpc.utah.edu/common/home/strong-group4/savanna/wrf_facet_data/";
poly_coef       = 1:3;
grad            = "opg";
years           = 2001:2015;
poly_fit        = "fitnlm";
grad_dir        = "/uufs/chpc.utah.edu/common/home/strong-group4/savanna/himat_opg_poly/" + radius + "km/";
cluster_number  = 5;
months          = [1:12];

% wrf file; where i get my elevation and lat lon data for plotting
fdir    = "/uufs/chpc.utah.edu/common/home/strong-group6/HIMAT/WRF_output/2000/output/";
fname   = "wrfout_d03_2000-01-01_00:00:00";

% location where I saved the OPG coefficients
corr_dir = "/uufs/chpc.utah.edu/common/home/strong-group4/savanna/himat_opg_poly/";

% precip data location
precip_dir  = "/uufs/chpc.utah.edu/common/home/strong-group4/savanna/precip_wrf/";

clust_dir   = "/uufs/chpc.utah.edu/common/home/strong-group4/savanna/cluster/clustStruct/";
% clust_dataset = "k-means_opg_coef_" + poly_coef_name + "_z-score_value_z-score_25km_fitnlm_nov19_bottomElevationZero.mat";
clust_dataset = "k-means_opg_coef_" + string(value) + "_value_" + string(value) + "_25km_fitnlm_nov19_bottomElevationZero_mean.mat";



%% Load clusters

load(clust_dir + clust_dataset);
clusters = cluster_data.clusters;
clusters = clusters(:, cluster_number);

facet_fname = "wrf_facets_" + radius;
load(facet_dir + facet_fname);


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

    for fi = 1:max(facets_labeled(:))
        count = 1;
        for mx = months
            polyVar_2D_monthlyx(fi, count) =  nanmean(polyVar_2D(fi, d_vec(:,2)==mx));
            count = count + 1;
        end
    end

    opg_opg_monthly = cat(2, opg_opg_monthly, polyVar_2D_monthlyx);
end

%% Fomrulate monthly average precipiation at each facet

avg_pr = [];

% Elevation
elev = ncread(fdir + fname, 'HGT', [1 1 1], [inf inf 1]); elev = elev';
elev = elev ./ 1000;

% loop to load in precip and coefficients
for yr = years
    disp(string(yr));
    
    % load in precip and dates
    load(precip_dir + "daily_precip_" + string(yr) + ".mat");
    numDays = size(t,1);
    % load coefficients
    load(corr_dir + radius + "km/daily_opg_" + string(yr) + "_" + radius + "km_fitnlm_nov19_bottomElevationZero.mat");
    
    avg_prx = NaN(max(facets_labeled(:)), numDays);

    % pull out daily values for the facets
    P1 = NaN(max(facets_labeled(:)), numDays);

    for fi = 1:max(facets_labeled(:))
        %disp(string(fi));
        use = facets_labeled == fi;
        delX = max(elev(use)) - min(elev(use)); % pull elevation change experienced on that facet
        eluse = elev(use) - min(elev(use)); % full elevation datapoints on the facet; mine make the lowest elevation 0 meters
    

        % loop through each day
        for dayx = 1:numDays

            if ~isnan(mode(opg_polyfit2(1,dayx,use)))
                prc = pr(dayx,use); % pull precip on that day on that facet
                [eluse_sort, idx] = sort(eluse); % sort elevation
                int_sum = trapz(eluse_sort, prc(idx))/delX; % calculate the facet average precip based on elevation change
            
                avg_prx(fi, dayx) = int_sum;
            end
        end

    end
    avg_pr = cat(2, avg_pr, avg_prx);

end


%% formulate monthly average based on cluster

clust_mean_pr = NaN(5, 12);

for cx = 1:5
    % pull facets at that cluster
    fi_pr = avg_pr(clusters==cx, :);

    fi_pr_mx = NaN(size(fi_pr, 1), 12);

    for mx = 1:12
        fi_pr_mx(:, mx) = mean(fi_pr(:, d_vec(:,2)==mx), 2, 'omitnan');

    end

    clust_mean_pr(cx, :) = mean(fi_pr_mx, 1, 'omitnan');

end

%% formulate clust logical 

clust_logical = [];
for cx = 1:5
    clust_logical = [clust_logical, clusters==cx];

end

clust_logical = logical(clust_logical);


%% plot the mean precipitation of each cluster and the median 
% precipitation total at the lowest elevation (p3) of each cluster

fig1 = figure('Position',[10 10 300 300]);

tt = tiledlayout(2,1);
tt.TileSpacing = 'compact';
tt.Padding = 'loose';

% Mean Precipiation
nexttile
for cx = 1:cluster_number

    plot(1:12, clust_mean_pr(cx, :), color=faceColor(cx), LineWidth=2);
    
    xlim([1,12]);
    xticks([1 2 3 4 5 6 7 8 9 10 11 12]);
    set(gca, 'XTickLabel', {'J' 'F' 'M' 'A' 'M' 'J' 'J' 'A' 'S' 'O' 'N' 'D'}, ...
        'fontsize', 7);
    ylim([0,7]);
    title("a) Precip. Average", 'fontsize', 7);
    ylabel("mm");
    hold on;
    grid on;
end

% Tile 3
nexttile
for cx = 1:cluster_number
    opg_median = nanmedian(opg_opg_monthly(clust_logical(:,cx), 25:36));

    disp(corrcoef(opg_median, clust_mean_pr(cx, :)));

    plot(1:12, opg_median, color=faceColor(cx), LineWidth=2);
    
    xlim([1,12]);
    xticks([1 2 3 4 5 6 7 8 9 10 11 12]);
    set(gca, 'XTickLabel', {'J' 'F' 'M' 'A' 'M' 'J' 'J' 'A' 'S' 'O' 'N' 'D'}, ...
        'fontsize', 7);
    ylim([0,5.5]);
    title("b) Precip. Total at the Facet's Lowest Elevation (P3)", 'fontsize', 7);
    ylabel("mm");
    hold on;
    grid on;
end

% create a legend
qw{1} = plot(nan, color=faceColor(1), LineWidth=3);
qw{2} = plot(nan, color=faceColor(2), LineWidth=3);
qw{3} = plot(nan, color=faceColor(3), LineWidth=3);
qw{4} = plot(nan, color=faceColor(4), LineWidth=3); 
qw{5} = plot(nan, color=faceColor(5), LineWidth=3); 
legend([qw{:}], {'Cluster 1','Cluster 2','Cluster 3', 'Cluster 4', 'Cluster 5'}, ...
    'location', 'southoutside', 'NumColumns', 3)

title(tt, "Median of Monthly OPG Coefficients", 'fontsize', 11);

% savefigname = ["/uufs/chpc.utah.edu/common/home/u1324060/himat_ms/pub_fig/" +...
%         grad + "_k-means_" + poly_name + "_for_" + cluster_number + "_median_coef.png"];
% export_fig(savefigname{:}, '-transparent', '-r500');


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
        rgb('red');];

    cmap = C(num,:);
    return
end




