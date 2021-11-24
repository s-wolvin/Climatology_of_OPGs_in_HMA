%% Savanna Wolvin
% created: Apr. 20th, 2021
% edited: Nov. 23rd, 2021

% SUMMARY
% pull the wrf output files and output daily varible values

% INPUT
% wrf_dir       - Directory holding the WRF output files
% save_dir      - Directory to save the WRF data as .mat files
% yr_initail    - Inital year of the WRF output files
% yr_final      - Final year of the WRF output files
% wrf_nest      - Define the WRF model nest to pull from
% wrf_var_name  - Define the WRF model variable to pull
% start_time    - Start and end time of 24-hr data to pull
% step_size     - Define number of timesteps within a 24-hr period

% OUTPUT
% WRF-NEST_WRF-VAR-NAME_YEAR_24h_05Z.mat - .mat file holding the daily WRF 
%                   model values and date vectors




%% clear workspace
clear; close all;




%% variable presets

% WRF output data location
wrf_dir = "/uufs/chpc.utah.edu/common/home/strong-group6/HIMAT/WRF_output/";

% Save folder
save_dir = "/uufs/chpc.utah.edu/common/home/strong-group4/savanna/var_wrf/";

% start and end year
yr_initial = 2001;
yr_final = 2015;

wrf_nest = "wrfout_d03";

wrf_var_name = "PSFC";

start_time = 5; % d03: 5, d02: 6
step_size = 24; % d03: 24, d02: 8

yr = yr_initial:yr_final;




%% create variables

wrf_fname = "2001/output/" + wrf_nest + "_2001-01-01_00:00:00";
wrf_file  = wrf_dir + wrf_fname;

LAT = ncread(wrf_file, 'XLAT',[1 1 1],[inf inf 1]); LAT = double(LAT);
LON = ncread(wrf_file, 'XLONG',[1 1 1],[inf inf 1]); LON = double(LON);

% creates number matrix for dates in WRF run
tref = (datenum([yr_initial 1 1 start_time 0 0]):datenum([yr_final+1 1 1 start_time 0 0]))';

% create a length X 402 X 525 NaN matrix
varx = NaN(length(tref),size(LAT,2),size(LON,1));




%% Loop through each year and pull daily averages
for yrx = yr

    file_path = wrf_dir + num2str(yrx) + "/output/";
    
    % pull all the output files in inner nested run
    di = dir(file_path + wrf_nest + "_" + yrx + "*");
    
    % pull each individual WRF output file under that year
    for i=1:length(di)
        idx = i + 1;
        fname = di(i).name;
        disp(fname);
        
        % skip loop b/c its the wrong file
        if strcmp(fname(end),'_')
            continue
        end
        
        % read the time values for one month of the year & convert to nums
        wtm = ncread(file_path + fname,'Times')';
        wtm = datenum(wtm);
        
        % if the wtm hour doesnt start at hour 0, skip loop
        f = find(hour(wtm)==start_time,1);
        if isempty(f)
            continue;
        end
        
        % pull timestamp for each day of 31 days
        wtm05Z = wtm(f:step_size:end);
        index_hr = (f:step_size:length(wtm))';
        dvec = datevec(wtm05Z);
        t_len = length(wtm);
        
        WRF_VAR = NaN(size(wtm05Z,1), size(LAT,2),size(LON,1));
        
        % daily loop
        for dayx = 1:length(index_hr)
            
            if (dayx == length(index_hr)) && (dvec(1,2) ~= 12)
                % pull current month
                first_range = (t_len+1) - index_hr(dayx);
                wvx_first = ncread(file_path + fname,wrf_var_name ,[1 1 index_hr(dayx)],[inf inf first_range],[1 1 1]);
                wvx_first = permute(wvx_first,[3 2 1]);
                
                % pull following month
                fname2 = di(idx).name;
                wtm2 = ncread(file_path + fname2,'Times')';
                f2 = find(hour(wtm2)==start_time,1); % pull indecies of the next file

                wvx_second = ncread(file_path + fname2,wrf_var_name ,[1 1 1],[inf inf (f2-1)],[1 1 1]);
                wvx_second = permute(wvx_second,[3 2 1]);
                
                % combine
                wvx = [wvx_first; wvx_second];
                wvx = squeeze(mean(wvx, 1));
                WRF_VAR(dayx,:,:) = wvx;
                
                
            elseif (dayx == length(index_hr)) && (dvec(1,2) == 12) && (dvec(1,1) ~= yr_final)
                % pull last month in year
                first_range = (t_len+1) - index_hr(dayx);
                wvx_first = ncread(file_path + fname,wrf_var_name,[1 1 index_hr(dayx)],[inf inf first_range],[1 1 1]);
                wvx_first = permute(wvx_first,[3 2 1]);
                
                % pull month of next year
                file_path2 = wrf_dir + num2str(yrx+1) + "/output/";
                di2 = dir(file_path2 + wrf_nest + "*");
                fname2 = di2(1).name;
                
                wtm2 = ncread(file_path2 + fname2,'Times')';
                f2 = find(hour(wtm2)==start_time,1); % pull indecies of the next file

                wvx_second = ncread(file_path2 + fname2,wrf_var_name,[1 1 1],[inf inf (f2-1)],[1 1 1]);
                wvx_second = permute(wvx_second,[3 2 1]);
                
                % combine
                wvx = [wvx_first; wvx_second];
                wvx = squeeze(mean(wvx, 1));
                WRF_VAR(dayx,:,:) = wvx;

            elseif (dayx == length(index_hr)) && (dvec(1,2) == 12) && (dvec(1,1) == yr_final)
                % last file of all years
                last_range = (t_len+1) - index_hr(dayx);
                wvx = ncread(file_path + fname,wrf_var_name,[1 1 index_hr(dayx)],[inf inf last_range],[1 1 1]);
                wvx = permute(wvx,[3 2 1]);
                
                wvx = squeeze(mean(wvx, 1));
                WRF_VAR(dayx,:,:) = wvx;
                
                
            else
                % regular loopvarx
                wvx = ncread(file_path + fname,wrf_var_name,[1 1 index_hr(dayx)],[inf inf step_size],[1 1 1]);
                % *** TIME X LAT X LON ***
                wvx = permute(wvx,[3 2 1]);
        
                % take average 
                wvx = squeeze(mean(wvx, 1));
                WRF_VAR(dayx,:,:) = wvx;
            end
        end
       
        use = ismember(tref,wtm);
        varx(use,:,:) = WRF_VAR;
    end
end




%% Loop through each year and save .mat files
varx = varx(1:end-1,:,:);
tref = tref(1:end-1);

for yrx = yr
    disp(yrx); 
    
    % of tref, which days are in the year we are saving?
    use = year(tref)==yrx; 
    
    % t will contain the datanum of the year we are saving 
    t = tref(use);
    
    % pr will contain the precip of the year we are saving
    wrf_var = varx(use,:,:);
   
    % save the precip and time datenum for that year
    save([save_dir + wrf_nest + "_" + wrf_var_name + "_" + num2str(yrx) + "_24h_05Z.mat"],'wrf_var','t', 'LAT', 'LON', '-v7.3');    

    % ????? all times and 10 x 10 lat-lon
    m = squeeze(WRF_VAR(:,10,10)); 
    datevec(t(m<0))
end



