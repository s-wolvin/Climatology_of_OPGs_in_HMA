## Savanna Wolvin
# Created: Feb. 26th, 2021
# Edited: Jul. 11th, 2022
    
# SUMMARY
# This program is used to pull 3D variables from the WRF Output files 
# and turn them into Mat Files. WRF output files for each year are in separate 
# files of that year number. This code can loop throgh each output folder year 
# and each monthly file within. 

# INPUT
# start_yr          - Start year of the whole simulation
# end_yr            - End year of the WRF simulation
# starting hour     - UTC hour to pull data
# filelocation      - Location of WRF files
# savelocation      - Location to save data
# threed_var        - 3D variable name 
#                       https://wrf-python.readthedocs.io/en/latest/diagnostics.html
# var_name          - What to name variable in MAT file

# OUTPUT
# [matfilename].mat - Output file of 3D variable for timeseries specified


#%% Imports
import numpy as np
import scipy.io as sio
import netCDF4 as nc
from wrf import getvar, ALL_TIMES, extract_times
import os
from datetime import timedelta, datetime
import pandas as pd
from numpy import arange


#%%
def main():
    #%% Variable Preset
    
    start_yr = 2001
    end_yr = 2015
    
    starting_hour = '05'
    
    #file location
    filelocation = "/uufs/chpc.utah.edu/common/home/strong-group6/HIMAT/WRF_output/"
    
    t_step   = 24 # D01: 2, D02: 8, D03: 24
        
    lat_size = 402 # D01: 114, D02: 219, D03: 402
    lon_size = 525 # D01: 145, D02: 270, D03: 525
    
    #directory path to where the figure will be saved
    savelocation = "/uufs/chpc.utah.edu/common/home/strong-group4/savanna/var_wrf/"
    
    # 3D vars: PWAT: 'pw', 2-m RH: 'rh2', 2-m TD: 'td2', 
    threed_var = 'pw'
    
    # what to name the variable in the new saved file
    var_name = 'PWAT'
    
    
    #%% load data
    
    years = np.arange(start_yr,(end_yr+1))
    
    for yr in years:
        file_path = filelocation + str(yr) + "/output/"
    
        files = os.listdir(file_path)
        d03Files = sorted([i for i in files if "wrfout_d03_" in i])
        d03Files = [i for i in d03Files if "." not in i] # remove any bad files
        
        tref = np.arange(datetime(yr, 1, 1, int(starting_hour)), datetime(yr+1, 1, 1, int(starting_hour)), timedelta(days=1), dtype='object')
        tref = pd.to_datetime(tref)
        
        # empty variable to hold a years worth of data
        VAR = np.zeros((np.shape(tref)[0], lat_size, lon_size))
        
        timecount = 0
        timecountend = 0
        
        for filex in range(len(d03Files)):
            print(d03Files[filex])
            
            #open netCDF file (ncread)
            ncfile = nc.Dataset(file_path + d03Files[filex], 'r')
            
            #turns array into normal array (masked ones suck)
            ncfile.set_always_mask(False) 
    
            # read in time value
            t_hour, t_year, t_mont, t, t_len = pullYearHour(ncfile)
                
            # pull index values for the starting hour
            index_05Z = [i for i in range(len(t_hour)) if (t_hour[i] == starting_hour) and (t_year[i] == str(yr))]
            
            # pull timestamps for each day
            t = t[index_05Z]
            
            #open the precipitable water variable
            #VARX = getvar(ncfile,threed_var,timeidx=index_05Z, meta=False)
            LAT = getvar(ncfile,'XLAT',timeidx=0, meta=False)
            LON = getvar(ncfile,'XLONG',timeidx=0, meta=False)
            
            rawdata = np.zeros((np.shape(index_05Z)[0], np.shape(LAT)[0], np.shape(LAT)[1]))
            
            for dayx in range(len(index_05Z)):
                timex = np.arange(index_05Z[dayx], (index_05Z[dayx]+t_step))
                
                if (dayx == (len(index_05Z)-1)) & (t_year[len(t_year)-1] == str(yr)) & ((t_year[len(t_year)-1] != str(end_yr)) or (t_mont[0] != '12')):
                    # pull following  month of data
                    first_range = np.arange(index_05Z[dayx], t_len)
                    first_halfd = getvar(ncfile, threed_var, timeidx=first_range, meta=False)
                    
                    second_halfd = secondFile(file_path, d03Files[filex+1], starting_hour, yr, threed_var)
                    
                    var  = np.concatenate((first_halfd, second_halfd), axis = 0)
                    
                    rawdata[dayx,:,:]= np.mean(var, axis = 0)
                    
                    
                elif (dayx == (len(index_05Z)-1)) & (t_year[len(t_year)-1] != str(yr)) & (t_year[len(t_year)-1] != str(end_yr+1)):
                    # Pull following year and month of data
                    first_range = np.arange(index_05Z[dayx], t_len)
                    first_halfd = getvar(ncfile, threed_var, timeidx=first_range, meta=False)
                    
                    # pull from the folder for the following year
                    file_path_next = filelocation + str(yr+1) + "/output/"
    
                    files_next = os.listdir(file_path_next)
                    d03Files_next = sorted([i for i in files_next if "wrfout_d03_" in i])
                    d03Files_next = [i for i in d03Files_next if "." not in i] # remove any bad files
                    
                    second_halfd = secondFile(files_next, d03Files_next[0], starting_hour, yr, threed_var)
                    
                    var  = np.concatenate((first_halfd, second_halfd), axis = 0)
                    
                    rawdata[dayx,:,:]= np.mean(var, axis = 0)
                    
                elif (dayx == (len(index_05Z)-1)) & (t_mont[0] == '12') & (t_year[len(t_year)-1] == str(end_yr)):
                    # pull the last few hour of data
                    first_range = np.arange(index_05Z[dayx], t_len)
                    first_halfd = getvar(ncfile, threed_var, timeidx=first_range, meta=False)
                    rawdata[dayx,:,:]= np.mean(first_halfd, axis = 0)
                    
                    
                else:
                    #open the vertical coordinate file            
                    rawdata[dayx,:,:]= np.mean(getvar(ncfile, threed_var, timeidx=timex, meta=False), axis = 0)
            
    
            timecountend += len(t)
            
            VAR[timecount:timecountend, :, :] = rawdata[:, :, :]
            
            timecount += len(t)
    
            #close netcdf file when finished
            ncfile.close()
    
    
        matfilename = 'wrfout_d03_' + threed_var + '_' + str(yr)
        
        #create the data dictonary
        data = {var_name:VAR,'LAT':LAT,'LON':LON}
    
        #save data dictonary to .mat file
        sio.savemat(savelocation + matfilename + '.mat', data)



def pullYearHour(NCFILE):    
    t_object = extract_times(NCFILE, timeidx=ALL_TIMES)
    t_datetime = pd.to_datetime(t_object)
    t_len = len(t_datetime)
    
    t_hour_a = t_datetime.strftime('%H')
    t_mont_a = t_datetime.strftime('%m')
    t_year_a = t_datetime.strftime('%Y')
    
    t_hour_a = t_hour_a.to_list()
    t_mont_a = t_mont_a.to_list()
    t_year_a = t_year_a.to_list()
    
    return t_hour_a, t_year_a, t_mont_a, t_datetime, t_len;


def pullIndex05(starting_hour, yr, t_hour, t_year):
    idx = [i for i in range(len(t_hour)) if (t_hour[i] == starting_hour) and (t_year[i] == str(yr))]
    
    return idx;


def secondFile(path, file, start_hour, year, threed_var):
    ncfile2 = nc.Dataset(path + file, 'r')
    t_hour2, t_year2, _, _, _ = pullYearHour(ncfile2)
    
    index_05Z2 = pullIndex05(start_hour, year, t_hour2, t_year2)
    
    second_range = arange(0, index_05Z2[0])
    second_var = getvar(ncfile2, threed_var, timeidx=second_range, meta=False)
    
    if (np.ndim(second_var) != 4):
        second_var  = np.expand_dims(second_var , axis = 0)
    
    ncfile2.close()

    return second_var




main()



