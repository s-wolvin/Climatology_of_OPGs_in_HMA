%% Savanna Wolvin
% Created: Aug. 30th, 2022
% Edited: Oct. 6th, 2022

% SUMMARY
% Conduct k-means clustering of OPG coefficients and save a data file
% containing the labeled clusters, centroids, and distances


%% Add Paths

addpath('/uufs/chpc.utah.edu/common/home/strong-group4/savanna/m_map');
addpath('/uufs/chpc.utah.edu/common/home/strong-group4/savanna/countries_shapefile');
addpath('/uufs/chpc.utah.edu/common/home/u1324060/exp_fig');



%% User-defined values

max_cluster     = 10;

years           = 2001:2015;

radius          = "25";
poly_fit        = "fitnlm";
grad            = "opg";
grad_dir        = "/uufs/chpc.utah.edu/common/home/strong-group4/savanna/himat_opg_poly/" + radius + "km/";
poly_coef       = 1:3;
poly_name       = "z-score_timeseries";

% raw, z-score_timeseries, z-score_facets, z-score_all, normalizeX-Axis,
% median_timeseries
% value           = "z-score_timeseries";
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


%% Colorbar

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
    rgb('midnightblue');...
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
    rgb('midnightblue');...
    rgb('gold');...
    rgb('black');];



%% Cluster

distance = NaN(max_cluster+1, max_cluster);
centroids = cell(max_cluster, 1);
clusters = [];

fig = figure;
% fig = figure('Visible','off');
rng(0);

% for seed = 1:100
%     rng(seed);

for cnum = 1:max_cluster
    % cnum = 5;
    disp(string(cnum));
    [idx, centroid, sumD] = kmeans(opg_opg_monthly, cnum);

    % decending order
    [idx,cluster_counts, sumD, centroid] = relabel_clusters(idx, sumD, centroid);

    % arrays to save
    distance(1:cnum, cnum) = sumD;
    centroids{cnum, 1} = centroid;
    clusters = cat(2, clusters, idx);

    % rework the data
    cluster_tempPlot = nan(size(facets_labeled));
    
    for fi = 1:length(idx)
        use = (facets_labeled == fi);
        cluster_tempPlot(use) = idx(fi);
    end

    % create colorbar
    Cx = C(1:cnum,:);

    % start cluster plot
    m_proj('sinusoidal','lat',[23.2,38.2],'lon',[65.5,86.4],'rectbox','on',  0.5);
    hold on
    m_coast('patch',[.7 .7 .7],'edgecolor','none'); 
    m_grid('box','on','xaxislocation', 'top','FontSize',12,'linestyle',':');
    
    hold on
    
    color = m_pcolor(lon, lat, cluster_tempPlot);
    set(color, 'facealpha', 0.8);
    colormap(Cx)
    c_line = rgb('DarkSlateGray');
    M=m_shaperead('ne_50m_admin_0_countries'); 
    for k=1:length(M.ncst)
         m_line(M.ncst{k}(:,1),M.ncst{k}(:,2),'color', c_line, ...
             'linewidth', 2);
    end
    
    set(gca,'FontSize',12)
    set(gcf,'Renderer','zbuffer')
    %title(string(cnum) + " Clustered Regions")
    
    caxis([0.5,cnum+0.5]);
    colorbar('Ticks', [1:cnum], 'fontsize', 12);


    savefigname = [fig_folder + grad + "_k-means_map_for_" +...
        num2str(cnum) + "_coef_" + poly_name + ".png"];
    savefigname = ["/uufs/chpc.utah.edu/common/home/u1324060/himat_ms/pub_fig/" +...
        grad + "_k-means_map_for_" +...
        num2str(cnum) + "_coef_" + poly_name + "_median.png"];
    % export_fig(savefigname{:}, '-transparent', '-r500');

    disp(string(sumD));

    clf;
    % 
    sse = nansum(distance.*distance,1);
    plot(1:cnum, sse(1:cnum)', "-ob", 'MarkerFaceColor', 'b');
    xlim([0,(cnum+1)]);
    xlabel("Number of Clusters");
    ylabel("Sum of the Squared Distance");
    grid('on');
    ylim([0,max(sse)*1.1]);

    savefigname = [fig_folder + grad + "_k-means_sumDistance_coef_" +...
        poly_name + "_median.png"];
    % export_fig(savefigname{:}, '-transparent', '-r500');

    clf;
% 
%     scatter(1:cnum, distance(:,1:cnum), 15, 'filled', 'MarkerFaceColor', 'b');
%     xlim([0,(cnum+1)]);
%     set(gca, 'YScale', 'log');
%     xlabel("Number of Clusters");
%     ylabel("Sum Distance");
% 
%     savefigname = [fig_folder + grad + "_k-means_sumDistance_log_coef_" +...
%         poly_name + ".png"];
%     export_fig(savefigname{:}, '-transparent', '-r500');

    % clf;


end

cluster_data.clusters  = clusters;
cluster_data.centroids  = centroids;
cluster_data.distance  = distance;
cluster_data.poly_coef = poly_coef;
cluster_data.months    = months;
cluster_data.var_type  = value;

save_file1 = data_folder + "k-means_opg_coef_" + poly_name +...
    "_value_" + string(value) + "_" + radius +...
    "km_fitnlm_nov19_bottomElevationZero_median.mat";
% save([save_file1],'cluster_data');

