%% Savanna Wolvin
% Created: Jun. 14th, 2021
% Edited: Jul. 21st, 2022

% SUMMARY
% Pulls the precipitation data and the facet data, fits the elevation
% -precipitation gradient on each facet with a linear and quadratic model.
% This script saves the elevation-precipitation gradient coefficients of 
% the quadratic and linear model, as well as their p-values, r-squared
% values, and AICc values.

% INPUT
% radius_dist -     Radius or radii used for faceting the terrain from
%                   wrf_facets_MAIN.m
% years -           Years of precipitation data to loop through
% p_thres -         The average facet precipiation needed for it to be 
%                   defined as a preciaption event
% wrf_dir -         Directory of the WRF output for the elevation data
% wrf_fname -       File name of one WRF output file
% precip_dir -      Directory of the precipitation data
% precip_file -     File name of the precipitation data
% save_dir -        Directory to save to

% OUTPUT
% daily_opg_" + num2str(yr) + "_" + radius_dist + "km_fitnlm_nov19_bottomElevationZero.mat
% .mat file with the quadratic coefficients

% daily_opg_" + num2str(yr) + "_" + radius_dist + "km_fitnlm_stats_nov19_bottomElevationZero.mat"
% .mat file with the p-value for all three coefficients, the r-squared 
% value, and the AICc index of the quadratic fits

% daily_opg_" + num2str(yr) + "_" + radius_dist + "km_fitlm_nov19_bottomElevationZero.mat
% .mat file with the linear coefficients

% daily_opg_" + num2str(yr) + "_" + radius_dist + "km_fitlm_stats_nov19_bottomElevationZero.mat
% .mat  file with the r-squared value and the AICc index of the linear fits




%% Variable presets
% radius distance of the faceting algorithm
radius_dist = "25";

years = 2001:2015;

% Threshold for mean precipitation on the facet to be an event
p_thres = 0.1;

% wrf output file for elevation
wrf_dir   = "/uufs/chpc.utah.edu/common/home/strong-group6/HIMAT/WRF_output/2000/output/";
wrf_fname = "wrfout_d03_2000-01-01_00:00:00";
wrf_file  = wrf_dir + wrf_fname;

% facet data location
facet_dir   = "/uufs/chpc.utah.edu/common/home/strong-group4/savanna/wrf_facet_data/";

% precipitation data
precip_dir = "/uufs/chpc.utah.edu/common/home/strong-group4/savanna/precip_wrf/";
precip_file = precip_dir + "daily_precip_" + years + "_24h_05Z.mat";

% save files to this location
save_dir    = "/uufs/chpc.utah.edu/common/home/strong-group4/savanna/himat_opg_poly/";




%% Pull elevation data
disp("Load elevation...");

elev = ncread(wrf_file,'HGT',[1 1 1],[inf inf 1]);
elev = elev';




%% Determine which facets are flat ground to skip curve-fitting

th      = 5; % max angle for terrain to be flat. (degrees)
lati = ncread(wrf_file,'XLAT',[1 1 1],[inf inf 1]); lati = double(lati'); 
loni = ncread(wrf_file,'XLONG',[1 1 1],[inf inf 1]); loni = double(loni'); 

% smooth elevation by guassian, w/ STD of sigma
elev_smoothed = imgaussfilt(elev,40/3/4);  % 3-sigma is 30 km, so sigma = 30/3/4km

% calculate direction topo is facing
[dx, dy] = gradient(elev_smoothed);
dy = -dy;
dx = -dx; dy = -dy; % flip orientation

% determine what gridpoints are flat.
flat = abs(dx)<th & abs(dy)<th; % slope<=prctile(slope,60);




%% Loop through each desired facet sampling radii

count = 1;

for radii = radius_dist
    disp("Load facets at a sampling radius of " + radii + "...");

    % Load facets
    facet_fname = "wrf_facets_" + radii;
    load(facet_dir + facet_fname); facets = facets';

    num_of_facets = unique(facets_labeled(:));
    facet_loop = length(num_of_facets);
    % facets_labeled(321,517) = 1533;




    % Loop through each year
    for yr = years
        disp(string(yr));

        % Pull precipitation
        load(precip_file(count),'pr','t');
        pr(pr<0) = NaN; % set any precip value less than zero to NaN

        % Create empty matricies to hold model coefficients, and stats
        % Quadratic Variables
        opg_polyfit2    = NaN(3,size(pr,1),size(pr,2),size(pr,3)); 
        p_value_quad    = NaN(3,size(pr,1),max(num_of_facets));
        r_squared_quad  = NaN(size(pr,1),max(num_of_facets)); 
        aicc_quad       = NaN(size(pr,1),max(num_of_facets));
        
        % Linear Variables
        opg_polyfit1   = NaN(2,size(pr,1),size(pr,2),size(pr,3)); 
        r_squared_lin  = NaN(size(pr,1),max(num_of_facets));
        aicc_lin       = NaN(size(pr,1),max(num_of_facets));

        precip_loop = size(pr,1);




        % Loop through each day of that year
        for dy = 1:precip_loop
            pd = squeeze(pr(dy,:,:)); % pull precip that one day

            % Fit OPG for that day for each facet
            for fi = 1:facet_loop
                use = facets_labeled == num_of_facets(fi);
                isflat = sum(flat(use));

                if isflat > 0
                    continue; % Facet is flat and is skipped
                end
                

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % LOOP TO REMOVE FACETS WHICH RESIDED OVER CENTRAL INDIA %%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if max(lati(use)) < 29
                    if sum(fi == [15, 149 243, 361]) ~= 1
                        bw = bwboundaries(use); bw = bw{1,1}; lat_facet = lati(bw(:,1));

                        bw2 = bwboundaries(facets_labeled == 996); bw2 = bw2{1,1};
                        lat_border = bw2(1587:2387,:);

                        idx = (min(bw(:,2)) <= lat_border(:,2)) .* (max(bw(:,2)) >= lat_border(:,2));
                        border_lat = 0;
                        for i = 1:length(idx)
                           if idx(i) == 1
                               border_lat = border_lat + lati(lat_border(i, 1), lat_border(i, 2));
                           end
                        end

                        if max(lat_facet) < ((border_lat/sum(idx))-0.5)
                            continue;
                        end
                    end
                end
                

                % Pull precipitation and elevation for this facet
                puse = pd(use);  % precip at facet 1, day 1
                eluse = elev(use);

                % Remove grid points with NaN precip values
                if sum(isnan(puse)) > 0
                    nan_index = isnan(puse);
                    eluse = eluse(~nan_index);
                    puse = puse(~nan_index);               
                end          

                psum = sum(puse);
                pmean = mean(puse, 'omitnan');
                plength = length(puse);

                % Skip facet if it have less than 3 grid points or is below
                % precipitation threshold
                if plength<3 || pmean < p_thres
                    continue;
                end

                % Shift facet elevations so that the lowest value is zero
                eluse = eluse - min(eluse(:));
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%% QUADRATIC MODEL
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                modelfun = @(b,x) b(1)*x.^2 + b(2)*x + b(3);
                beta0 = [0,0,0];
                md1 = fitnlm((eluse./1000), puse, modelfun, beta0);

                % Save coefficients
                opg_polyfit2(1,dy,use) = md1.Coefficients{1,1};
                opg_polyfit2(2,dy,use) = md1.Coefficients{2,1};
                opg_polyfit2(3,dy,use) = md1.Coefficients{3,1};
                
                % Save stats
                p_value_quad(1,dy,fi) = md1.Coefficients{1,4};
                p_value_quad(2,dy,fi) = md1.Coefficients{2,4};
                p_value_quad(3,dy,fi) = md1.Coefficients{3,4};
                r_squared_quad(dy,fi) = md1.Rsquared.Ordinary;
                aicc_quad(dy,fi)      = md1.ModelCriterion.AICc;


                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%% LINEAR MODEL
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                modelfun = @(b,x) b(1)*x + b(2);
                beta0 = [0,0];
                md1 = fitnlm((eluse./1000), puse, modelfun, beta0);

                % Save coefficients
                opg_polyfit1(1,dy,use) = md1.Coefficients{1,1};
                opg_polyfit1(2,dy,use) = md1.Coefficients{2,1};

                % Save stats
                r_squared_lin(dy,fi) = md1.Rsquared.Ordinary;
                aicc_lin(dy,fi)      = md1.ModelCriterion.AICc;
                



            end
        end




        disp("Saving...");

        save_file1 = save_dir + radius_dist + "km/daily_opg_" + num2str(yr) + "_" + radius_dist + "km_fitnlm_nov19_bottomElevationZero.mat";
        save(save_file1,'opg_polyfit2','t');

        save_file2 = save_dir + radius_dist + "km/daily_opg_" + num2str(yr) + "_" + radius_dist + "km_fitnlm_stats_nov19_bottomElevationZero.mat";
        save(save_file2,'p_value_quad','r_squared_quad','aicc_quad','t');
        
        ave_file3 = save_dir + radius_dist + "km/daily_opg_" + num2str(yr) + "_" + radius_dist + "km_fitlm_nov19_bottomElevationZero.mat";
        save(save_file3,'opg_polyfit1','t');

        save_file4 = save_dir + radius_dist + "km/daily_opg_" + num2str(yr) + "_" + radius_dist + "km_fitlm_stats_nov19_bottomElevationZero.mat";
        save(save_file4,'r_squared_lin','aicc_lin','t');

        count = count + 1;

     end
end 



