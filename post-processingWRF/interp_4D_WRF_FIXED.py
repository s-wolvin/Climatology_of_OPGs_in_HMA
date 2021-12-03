## Savanna Wolvin
# Created: Feb. 25th, 2021
# Edited: Dec. 3rd, 2021
    
# SUMMARY
# This program is used to pull 4D variables from the WRF Output files 
# and turn them into Mat Files

# INPUT
# start_yr      - First year of the WRF output data to intepolate 
# end_yr        - Last year of the WRF output data to interpolate
# interp_method - What kind of levels do you want to interpolate to, such as a
#                   pressure level, or height level
# interp_levels - Which levels do you want to interpolate. I would suggest a 
#                   maximum of 4 levels due to the size of the file
# starting_hour - Since the interpolated data is saved as 24-hour daily 
#                   averages, therefore you can choose what hour in the 
#                   24-hour period the daily averages start
# fourd_var     - Which 4D variable to interpolate. only one at a time
# wrf_file      - Indicate when nest of the wrf model to use
# filelocation  - What is the file location? the WRF output files are expected 
#                   to be in yearly folders like $year/output/$wrf_files on 
#                   line 69, so only put file path up to the year
# savelocation  - Where to save the output mat files

# OUTPUT
# wrf_$fourd_var_700_600_500_400mb_$year.mat - mat file which holds the 
#                   interpolated atmospheric variables in the following order:
#                   (Levels, Days, Lat, Long)


#%% Imports
import numpy as np
import scipy.io as sio
import netCDF4 as nc
from wrf import getvar, interplevel, ALL_TIMES, extract_times
import os
from datetime import timedelta, datetime
import pandas as pd
from numpy import arange


#%% MAIN
def main():
    #%% Variable Preset
    
    start_yr = 2001
    end_yr = 2015
    
    #Interplation method
    #p = pressure, z = height, th = theta
    interp_method = 'p'
    interp_levels = [90000, 80000, 70000, 60000, 50000, 40000, 30000, 20000]
    
    # first_halfd = getvar(nc.Dataset(ilelocation + "2001/output/" + d03Files[filex], 'r'), 'ter', timeidx=0, meta=False)
    
    starting_hour = '06'
    
    # fourd_variables_to_interpolate = ['tk','td','ua','va','wa']
    fourd_var = 'geopt'
    
    wrf_file = "wrfout_d02_"
    
    lat_size = 219 # D01: 114, D02: 219, D03: 402
    lon_size = 270 # D01: 145, D02: 270, D03: 525
    
    #file location
    filelocation = "/uufs/chpc.utah.edu/common/home/strong-group6/HIMAT/WRF_output/"
    
    #directory path to where the figure will be saved
    savelocation = "/uufs/chpc.utah.edu/common/home/strong-group4/savanna/var_wrf/"
    
    
    #%% Load data
    
    years = np.arange(start_yr,end_yr+1)
    
    for yr in years:
        file_path = filelocation + str(yr) + "/output/"
    
        # load list of output files for this year
        files = os.listdir(file_path)
        d03Files = sorted([i for i in files if wrf_file in i])
        d03Files = [i for i in d03Files if "." not in i] # remove any bad files
        
        # create reference array of datetimes
        tref = np.arange(datetime(yr, 1, 1, int(starting_hour)), datetime(yr+1, 1, 1, int(starting_hour)), timedelta(days=1), dtype='object')
        tref = pd.to_datetime(tref)
        
        interpdata = np.zeros((np.shape(interp_levels)[0], np.shape(tref)[0], lat_size,lon_size))
        
        del tref
        
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
                   
            # index values of the starting hour
            index_05Z = pullIndex05(starting_hour, yr, t_hour, t_year)
            
            # pull timestamps for each day
            t = t[index_05Z]
            
            #open the lat lon variables
            LAT = getvar(ncfile,'XLAT',timeidx=0, meta=False)
            LON = getvar(ncfile,'XLONG',timeidx=0, meta=False)
            
            # create empty files for the vertical coordinates and atmospheric variables
            vertical_coord = np.zeros((np.shape(index_05Z)[0], 44 ,  np.shape(LAT)[0], np.shape(LAT)[1]))
            rawdata = np.zeros((np.shape(index_05Z)[0], 44 ,  np.shape(LAT)[0], np.shape(LAT)[1]))
            
            # loop through each day
            for dayx in range(len(index_05Z)):
                timex = np.arange(index_05Z[dayx], (index_05Z[dayx]+8))
                
                if (dayx == (len(index_05Z)-1)) & (t_mont[0] != '12'): 
                # pull following  month of data
                    first_range = np.arange(index_05Z[dayx], t_len)
                    first_halfv = getvar(ncfile,interp_method, timeidx=first_range, meta=False)
                    first_halfd = getvar(ncfile, fourd_var, timeidx=first_range, meta=False)
                    
                    second_halfv, second_halfd = secondFile(file_path, d03Files[filex+1], starting_hour, yr, interp_method, fourd_var)
                    
                    vert = np.concatenate((first_halfv, second_halfv), axis = 0)
                    var  = np.concatenate((first_halfd, second_halfd), axis = 0)
                    
                    del first_range, first_halfv, first_halfd, second_halfv, second_halfd
                    
                    vertical_coord[dayx,:,:,:] = np.mean(vert, axis = 0)
                    rawdata[dayx,:,:,:]= np.mean(var, axis = 0)
                    
                    del vert, var
                    
                elif (dayx == (len(index_05Z)-1)) & (t_mont[0] == '12') & (t_year[0] != str(end_yr)): 
                # Pull following year and month of data
                    first_range = np.arange(index_05Z[dayx], t_len)
                    first_halfv = getvar(ncfile,interp_method, timeidx=first_range, meta=False)
                    first_halfd = getvar(ncfile, fourd_var, timeidx=first_range, meta=False)
                    
                    file_path_next = filelocation + str(yr+1) + "/output/"
    
                    files_next = os.listdir(file_path_next)
                    d03Files_next = sorted([i for i in files_next if wrf_file in i])
                    d03Files_next = [i for i in d03Files_next if "." not in i] # remove any bad files
                    
                    second_halfv, second_halfd = secondFile(file_path_next, d03Files_next[0], starting_hour, yr+1, interp_method, fourd_var)
                    
                    vert = np.concatenate((first_halfv, second_halfv), axis = 0)
                    var  = np.concatenate((first_halfd, second_halfd), axis = 0)
                    
                    del first_range, first_halfv, first_halfd, second_halfv, second_halfd, files_next, d03Files_next
                    
                    vertical_coord[dayx,:,:,:] = np.mean( vert, axis = 0)
                    rawdata[dayx,:,:,:]= np.mean(var, axis = 0)
                    
                elif (dayx == (len(index_05Z)-1)) & (t_mont[0] == '12') & (t_year[0] == str(end_yr)):
                # pull the last few hour of data
                    first_range = np.arange(index_05Z[dayx], t_len)
                    first_halfv = getvar(ncfile,interp_method, timeidx=first_range, meta=False)
                    first_halfd = getvar(ncfile, fourd_var, timeidx=first_range, meta=False)
                    vertical_coord[dayx,:,:,:] = np.mean(first_halfv, axis = 0)
                    rawdata[dayx,:,:,:]= np.mean(first_halfd, axis = 0)
                    
                    del first_range, first_halfv, first_halfd
                    
                else:
                    #open the vertical coordinate file
                    vertical_coord[dayx,:,:,:] =  np.mean(getvar(ncfile,interp_method, timeidx=timex, meta=False), axis = 0)
                
                    rawdata[dayx,:,:,:]= np.mean(getvar(ncfile, fourd_var, timeidx=timex, meta=False), axis = 0)
            
            
            #%%
            count = 0
            timecountend += len(t)
            for level in interp_levels:
                #interpolate the data
                interpdata[count, timecount:timecountend, :, :] = interplevel(rawdata[:,:,:,:],vertical_coord,level, meta=False)
                count += 1
                # print(level)
                
            timecount += len(t)
            
            #close netcdf file when finished
            ncfile.close()
            
        matfilename = wrf_file + fourd_var + '_700_600_500_400mb_' + str(yr)
        
        #create the data dictonary
        data = {'LAT':LAT,'LON':LON}
    
        data.update({fourd_var[0]:interpdata[:,:,:,:]})
    
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


def secondFile(path, file, start_hour, year, interp_method, fourd_var):
    ncfile2 = nc.Dataset(path + file, 'r')
    t_hour2, t_year2, _, _, _ = pullYearHour(ncfile2)
    
    index_05Z2 = pullIndex05(start_hour, year, t_hour2, t_year2)
    
    second_range = arange(0, index_05Z2[0])
    second_vert = getvar(ncfile2,interp_method, timeidx=second_range, meta=False)
    second_var = getvar(ncfile2, fourd_var, timeidx=second_range, meta=False)
    
    if (np.ndim(second_var) != 4) & (np.ndim(second_vert) != 4):
        second_var  = np.expand_dims(second_var , axis = 0)
        second_vert = np.expand_dims(second_vert, axis = 0)
    
    ncfile2.close()

    return second_vert, second_var


main()









