%% Savanna Wolvin
% Created: Nov. 8th, 2020
% Edited: Jul. 21st, 2022

% SUMMARY
% Plot the elevation-precipitation relationship from user-defined facets
% and dates along with a histrogram illustrating the distribution of
% grid points with elevation.

% INPUT
% years -           Years of data to plot
% fnum -            Desired facet number
% radius -          Desired sampling radius of facet data
% wrf_dir -         Directory of the WRF output data
% wrf_fnmae -       WRF output file for elevation data
% facet_dir -       Directory of the facet data
% facet_fname -     Output files from wrf_facets_MAIN.m
% save_dir -        Directory to save to

% OUTPUT
% opg_' + num2str(fi) + '_' + datex + '_areaVSelev.png
% outputs transparent .png file with a scatter plot of precipitation VS 
% elevation, the quadratic curve-fit, and the average precipiaiton on the 
% facet found with respect to area and elevation. 




%% export figure path
addpath('/uufs/chpc.utah.edu/common/home/u1324060/exp_fig/');




%% variable pre-set


years = 2001:2015;

fnum = 939;

radius = "25";

% wrf data for elevation
wrf_dir   = "/uufs/chpc.utah.edu/common/home/strong-group6/HIMAT/WRF_output/2000/output/";
wrf_fname = "wrfout_d03_2000-01-01_00:00:00";

% facet data
facet_dir   = "/uufs/chpc.utah.edu/common/home/strong-group4/savanna/wrf_facet_data/";
facet_fname = "wrf_facets_" + radius;

% precip data
precip_dir   = "/uufs/chpc.utah.edu/common/home/strong-group4/savanna/precip_wrf/";

% save folder
save_dir = "/uufs/chpc.utah.edu/common/home/strong-group4/savanna/himat_opg_poly/";




%% load data
disp("Load data...");

% Load Facet data
load(facet_dir + facet_fname); 
facets = facets';
num_of_facets = unique(facets_labeled(:));
max_label = max(facets_labeled(:));

% Load Elevation data
elev = ncread(wrf_dir + wrf_fname,'HGT',[1 1 1],[inf inf 1]);
elev = elev';




%% Plot difference in averages


for yr = years
    disp(string(yr));

    % Load precip data
    precip_fname = "daily_precip_" + string(yr) + "_24h_05Z.mat";
    load(precip_dir + precip_fname);

    fig = figure('Position', [100,100,600,700]);
    



    for dayx = 1:size(pr,1)
        datex = datestr(t(dayx),29);
        
        for fi = fnum
            disp(string(fi));

            % Pull facet location
            use = (facets_labeled == fi);
            x = elev(use);
            
            % Skip small facets / facets with no precipitation
            if length(x) < 3 || mean(pr(dayx,use), 'omitnan') < 0.1 
                continue;
            end
            
            tt = tiledlayout(3, 2); 
            tt.TileSpacing = 'tight';
            tt.Title.FontWeight = 'bold';
            tt.Padding = 'loose';
            
            % Process facet elevations
            eluse = elev(use);
            eluse = eluse - min(eluse(:));
            elev_diff = (max(elev(use)) - min(elev(use))); 
            
            %%%%%%%% PRECIPITATION-ELEVATION PLOT %%%%%%%%%%%%%%%%%%%%%%%%%
            nexttile([2,2]);
            scatter(elev(use), pr(dayx, use), 6, [0.1215 0.2509 1], 'filled'); hold on;
            xlim([min(elev(use)), max(elev(use))]);
            mid = (max(pr(dayx, use)) + min(pr(dayx, use))) / 2;
            
            % AVERAGE BY AREA
            gpMean = mean(pr(dayx, use), 'omitnan');
            pp3 = plot([min(elev(use)), max(elev(use))], [gpMean,gpMean],...
                'LineWidth', 3);
            pp3.Color(4) = 0.9;
            pp3.Color = [1,0,0];
            text(max(elev(use)), gpMean-3, {string(round(gpMean,1)) + " mm"},...
                'fontsize', 10, 'color', 'r', 'HorizontalAlignment', 'right');
            
            % AVERAGE BY ELEVATION
            tmp = pr(dayx, use);
            [eluse_sort, idx] = sort(eluse);
            delX = max(elev(use)) - min(elev(use));
            tpMean = trapz(eluse_sort, tmp(idx))/delX;
            pp2 = plot([min(elev(use)), max(elev(use))], [tpMean,tpMean],...
                'LineWidth', 3);
            pp2.Color(4) = 0.9;
            pp2.Color = [0,0,0];
            text(max(elev(use)), tpMean-3, {string(round(tpMean,1)) + " mm"},...
                'fontsize', 10, 'HorizontalAlignment', 'right');
            
            % PLOT LABELS
            legend(["WRF Precip.","Area Avg.", "Elevation Avg."],...
                'location', 'northwest', 'fontsize', 11);
            title({"Precipitation on Facet " + num2str(fi); datestr(t(dayx), 22)});
            ylabel("24-Hour Precipitation Total (mm)", 'fontsize', 12);
            xticklabels([]);
            grid on;
            set(gca, 'fontsize', 12);
            
            %%%%%%%% ELEVATION HISTOGRAM PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            nexttile([1,2]);
            histogram(elev(use), min(elev(use)):(elev_diff/20):max(elev(use)),...
                'FaceColor',[0.1215 0.2509 1]);
            xlim([min(elev(use)), max(elev(use))]);
            ylabel("Grid Point Count");
            grid on;
            
            % TOTAL PLOT LABELS, ETC
            set(gca, 'fontsize', 12);
            xlabel(tt, "Elevation (m)", 'fontsize', 12);        
            
            hold off;

            save_file = save_dir + fi + '/opg_' + num2str(fi) + '_' + datex + '_areaVSelev.png';
            export_fig(fig, save_file{:}, '-transparent', '-r300');
    
            clf;
        end
    end

end



