%% Savanna Wolvin
% created: Mar. 18th, 2021
% edited: Nov. 23rd, 2021

% SUMMARY
% pull the wrf output files and output single hour values during a 
% specified hour

% INPUT
% wrf_dir       - Directory holding the WRF output files
% save_dir      - Directory to save the WRF data as .mat files
% yr_initail    - Inital year of the WRF output files
% yr_final      - Final year of the WRF output files
% hourx         - Single hour to pull data from
% hourIsMorningAfter - True to remove the first value instead of the last
%                       value

% OUTPUT
% hourly_WRF-VAR-NAME_YEAR_HOURZ.mat - .mat file holding a specified hour 
%                                          of a specified WRF Model
%                                          variable and date vectors




%% clear workspace
clear; close all;




%% variable presets

% WRF output data location
wrf_dir = "/uufs/chpc.utah.edu/common/home/strong-group6/HIMAT/WRF_output/";

% Save folder
save_dir = "/uufs/chpc.utah.edu/common/home/strong-group4/savanna/precip_wrf/";

% start and end year
yr_initial = 2001;
yr_final = 2015;

yr = yr_initial:yr_final;

hourx = 1; % max 8:10, min 0:1, 23

% cut off first hour??
hourIsMorningAfter = 1;

wrf_var_name = 'T2';




%% create variables

% creates number matrix for dates in WRF run
tref = (datenum([yr_initial 1 1 hourx 0 0]):datenum([yr_final+1 1 1 hourx 0 0]))';

% create a length X 402 X 525 NaN matrix
wrf_varX = NaN(length(tref),402,525);




%% loop through each year and pull single hour
for yrx = yr

    file_path = wrf_dir + num2str(yrx) + "/output/";
    
    % pull all the output files in inner nested run
    di = dir(file_path + "wrfout_d03*");
    
    % pull each individual WRF output file under that year
    for i=1:length(di)
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
        f = find(hour(wtm)==hourx, 1);
        if isempty(f)
            continue;
        end
        
        % pull timestamp for each day of 31 days
        wtm = wtm(f:24:end);
        
        % pull 2m temp
        varx = ncread(file_path + fname,wrf_var_name,[1 1 f],[inf inf inf],[1 1 24]);
        
        % *** TIME X LAT X LON ***
        varx = permute(varx,[3 2 1]);
        
        % check if wrf time is NOT memeber of wanted time, 0 for is member
        k = ~ismember(wtm,tref);
        % if timestamp isnt in tref, then emtpy that location
        wtm(k) = [];
        % if timestamp isnt in tref, then empty that precip location
        varx(k,:,:) = [];
        % find location of timestamp an put into precip matrix
        use = ismember(tref,wtm);
        wrf_varX(use,:,:) = varx;
        
 
    end
end




%% Save WRF Variable

if hourIsMorningAfter
    wrf_varX = wrf_varX(2:end,:,:);
else
    wrf_varX = wrf_varX(1:end-1,:,:);
end

tref = tref(1:end-1);

for yrx = yr
    disp(yrx); 
    
    % of tref, which days are in the year we are saving?
    use = year(tref)==yrx; 
    
    % t will contain the datanum of the year we are saving 
    t = tref(use);
    
    % pr will contain the precip of the year we are saving
    wrf_var = wrf_varX(use,:,:);
   
    % save the precip and time datenum for that year
    save([save_dir + "hourly_" + string(wrf_var_name) + "_" + num2str(yrx) + "_" + hourx + "Z.mat"],'wrf_var','t');
    
    % ????? all times and 10 x 10 lat-lon
    m = squeeze(wrf_var(:,10,10)); 
    datevec(t(m<0))
end



