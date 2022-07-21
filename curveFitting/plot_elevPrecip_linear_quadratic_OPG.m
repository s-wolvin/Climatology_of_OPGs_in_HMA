%% Savanna Wolvin
% Created: Nov. 8th, 2020
% Edited: Jul. 21st, 2022

% SUMMARY
% This code plots the elevation VS precipitation with the gaussian fit 

% INPUT
% year -        year of data to plot
% radius -      desired sampling radius of facet data
% fnum -        desired facet number
% start_dy -    the inital day of the year to start plotting
% wrf_file -    wrf output file for elevation data
% facet_file -  output files from wrf_facets_MAIN.m
% precip_file - output files from pull_himat_precip_MAIN.m
% opg_file -    output files from calc_himat_opg_gauss.m
% save_dir -    directory to save the figures

% OUTPUT
% opg_(fnum)_(year)-(month)-(day).png - outputs transparent .png file with 
%                                       a scatter plot of precipitation VS 
%                                       elevation and a line plot of the
%                                       Gaussian OPG on top. Plotted for an
%                                       entire year for a single facet



%% export figure path
addpath('/uufs/chpc.utah.edu/common/home/u1324060/exp_fig/');




%% How to plot the figure?

% 1 - Show AICC
% 2 - Show Coefficients
% 3 - Show only Facet # and Date
title_options = 2;

% 1 - Show Quadratic
% 2 - Show Linear
% 3 - Show Quadratic and Linear
% 4 - Show Only the WRF Precipitation
plot_options = 1;

% 1 - Visisble plot
% 2 - Invisible plot
show_options = 1;




%% variable pre-set


year = "2001";

fnum = 939;

start_dy = 1;
end_dy = 365;
radius = "25";

% wrf data for elevation
wrf_dir   = "/uufs/chpc.utah.edu/common/home/strong-group6/HIMAT/WRF_output/2000/output/";
wrf_fname = "wrfout_d03_2000-01-01_00:00:00";

% facet data
facet_dir   = "/uufs/chpc.utah.edu/common/home/strong-group4/savanna/wrf_facet_data/";
facet_fname = "wrf_facets_" + radius;

% precip data
precip_dir   = "/uufs/chpc.utah.edu/common/home/strong-group4/savanna/precip_wrf/";
precip_fname = "daily_precip_" + year + "_24h_05Z.mat";

% quadratic and linear opg data
opg_dir          = "/uufs/chpc.utah.edu/common/home/strong-group4/savanna/himat_opg_poly/" + radius + "km/";
opg_fname_poly   = "daily_opg_" + year + "_" + radius + "km_fitnlm_nov19_bottomElevationZero.mat";
opg_fname_lin    = "daily_opg_" + year + "_" + radius + "km_fitlm_nov19_bottomElevationZero.mat";

% stats file
opg_stat_lin   = "daily_opg_" + year + "_" + radius + "km_fitlm_stats_nov19_bottomElevationZero";
opg_stat_qua   = "daily_opg_" + year + "_" + radius + "km_fitnlm_stats_nov19_bottomElevationZero"; 

% save folder
save_dir = "/uufs/chpc.utah.edu/common/home/strong-group4/savanna/himat_opg_poly/";




%% load data
disp("Load data...");

% Facet data
load(facet_dir + facet_fname); 
facets = facets';
num_of_facets = unique(facets_labeled(:));
max_label = max(facets_labeled(:));

% Load coefficients
if plot_options == 1 || plot_options == 3
    load(opg_dir + opg_fname_poly);
    load(opg_dir + opg_stat_qua);
end

if plot_options == 2 || plot_options == 3
    load(opg_dir + opg_fname_lin);
    load(opg_dir + opg_stat_lin);
end

% Load Precip
load(precip_dir + precip_fname);

% Load Elevation
elev = ncread(wrf_dir + wrf_fname,'HGT',[1 1 1],[inf inf 1]);
elev = elev';




%% plot

if show_options == 1
    fig = figure('Position', [10, 10, 650, 500]);
elseif show_optons == 2
    fig = figure('Position', [10, 10, 650, 500], 'Visible', 'off');
end




for dayx = start_dy:end_dy 
    disp(string(dayx));
    datex = datestr(t(dayx),29);
    
    for fi = fnum
        disp(string(fi));

        % pull facet location and precip mean
        use = (facets_labeled == fi);
        y_precip = pr(dayx, use);
        p_mean = mean(y_precip, 'omitnan');

        if p_mean < 0.1; continue; end % skip little to no precip days
        
        % Elevation Data for this Facet
        x = elev(use);
        x1 = [min(x):max(x)];
        x1v2 = (x1-min(x))./1000;

        % Scatter Plot of WRF Precipitation
        scatter(x./1000, y_precip, 6,[0.1215 0.2509 1],'filled', 'DisplayName', 'WRF Precip.');
        hold on; 

        % Plot Linear Fit
        if plot_options == 2 || plot_options == 3
            ft = fit(double(x./1000), y_precip', 'poly1');
            pp2 = plot(ft, 'DisplayName', "Linear Fit"); 
            pp2.Color(4) = 0.9;
            pp2.LineWidth = 4;
            pp2.Color = [1, 0.4666, 0];
        end

        % Plot Quadratic Fit
        if plot_options == 1 || plot_options == 3
            p1 = mode(opg_polyfit2(1,dayx,use), 'all');
            p2 = mode(opg_polyfit2(2,dayx,use), 'all');
            p3 = mode(opg_polyfit2(3,dayx,use), 'all');
            pp1 = plot(x1/1000, p1*(x1v2.^2) + p2*(x1v2) + p3,...
                'color', 'k', 'DisplayName', 'Quadratic Fit');
            pp1.Color(4) = 0.9;
            pp1.LineWidth = 4;
        end

        % Select Title
        if title_options == 1
            title(["Facet " + num2str(fi) + ": " + datestr(t(dayx), 22) ;...
                "Quadratic AICc = " + string(round(aicc(dayx, fi))) + ", " + "Linear AICc = " + string(round(aicc_lin(dayx, fi)))]);
        elseif title_options == 2
            title({"Facet " + num2str(fi) + ": " + datestr(t(dayx), 22);...
                "P1 = " + round(p1,2) + ", P2 = " + round(p2,2) + ", P3 = " + round(p3,2)});
        elseif title_options == 3
            title("Facet " + num2str(fi) + ": " + datestr(t(dayx), 22));
        end

        % Create Labels
        xlabel("Elevation (km)");
        ylabel("24-Hour Precipitation Total (mm)");
        ylim([0,max(y_precip)+1]);
        xlim([min(elev(use)), max(elev(use))]/1000);
        legend('location', 'northwest');
        
        ax = gca; ax.FontSize = 14;
        
        hold off;

        save_file = save_dir + fi + '/opg_' + num2str(fi) + '_' + datex + '_05Z.png';
        export_fig(fig, save_file{:}, '-transparent', '-r300');

        clf;

    end
end



