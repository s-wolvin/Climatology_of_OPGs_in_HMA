## Savanna Wolvin
# Created: Mar. 25th, 2021
# Edited: Dec. 3rd, 2021
    
# SUMMARY
# This program is used to pull a meridional cross section from the WRF Output 
# files and turn them into Mat Files

# INPUT
# start_yr      - First year of the WRF output data to intepolate 
# end_yr        - Last year of the WRF output data to interpolate
# interp_method - What kind of levels do you want to interpolate to, such as a
#                   pressure level, or height level
# starting_hour - Since the interpolated data is saved as 24-hour daily 
#                   averages, therefore you can choose what hour in the 
#                   24-hour period the daily averages start
# fourd_var     - Which 4D variable to interpolate. only one at a time
# start_lat     - initiial latitude point
# end_lat       - final latitude point
# longitude     - longitude for the cross section
# fnum          - facet number for the save file
# filelocation  - What is the file location? the WRF output files are expected 
#                   to be in yearly folders like $year/output/$wrf_files on 
#                   line 69, so only put file path up to the year
# savelocation  - Where to save the output mat files

# OUTPUT
# wrf_$fnum_$fourd_var_mer_cs_$year.mat - mat file which holds the 
#                   interpolated atmospheric variables in the following order:
#                   (Levels, Days, Lat, Long)


#%% imports
import numpy as np
from numpy import arange, mean, concatenate
import scipy.io as sio
import netCDF4 as nc
#import netcdf-c/4.7.4 as nc
from wrf import getvar, vertcross, ALL_TIMES, extract_times, CoordPair, to_np
import os
from datetime import timedelta, datetime
import pandas as pd


#%% MAIN
def main():
    
    #%% Variable Preset
    
    start_yr, end_yr = 2001, 2015
    
    #Interplation method
    #p = pressure, z = height, th = theta
    interp_method = 'z'
    
    starting_hour = '05'
    
    # fourd_variables_to_interpolate = ['tk','td','ua','va','wa']
    fourd_var = 'va'
    
    # pick start and end, lat and lons
    start_lat = 23.6621
    end_lat = 37.5574
    longitude = 76.1108
    
    fnum = 870
    
    #file location
    filelocation = "/uufs/chpc.utah.edu/common/home/strong-group6/HIMAT/WRF_output/"
    
    #directory path to where the figure will be saved
    savelocation = "/uufs/chpc.utah.edu/common/home/strong-group4/savanna/var_wrf/"
    
    
    #%% Load data
    
    start_point = coordinates(start_lat, longitude)
    end_point = coordinates(end_lat, longitude)
    
    years = arange(start_yr,end_yr+1)
    
    for yr in years:        
        file_path = filelocation + str(yr) + "/output/"
    
        # load list of output files for this year
        files = os.listdir(file_path)
        d03Files = sorted([i for i in files if "wrfout_d03_" in i])
        d03Files = [i for i in d03Files if "." not in i] # remove any bad files
        
        # create reference array of datetimes
        tref = arange(datetime(yr, 1, 1, int(starting_hour)), datetime(yr+1, 1, 1, int(starting_hour)), timedelta(days=1), dtype='object')
        tref = pd.to_datetime(tref)
        
        Z = np.zeros((np.shape(tref)[0], 150))
        Data = np.zeros((np.shape(tref)[0], 150, 388))
        
        del tref
        
        timecount = 0
        timecountend = 0
        
        for filex in range(len(d03Files)):
            print(d03Files[filex])
            
            #open netCDF file (ncread)
            ncfile = nc.Dataset(file_path + d03Files[filex], 'r')
            
            #turns array into normal array (masked ones suck)
            # ncfile.set_always_mask(False) 
            
            # read in time value
            t_hour, t_year, t_mont, t, t_len = pullYearHour(ncfile)
                   
            # index values of the starting hour
            index_05Z = pullIndex05(starting_hour, yr, t_hour, t_year)
            
            # pull timestamps for each day
            t = t[index_05Z]
            
            Z_file = np.zeros((np.shape(index_05Z)[0], 150))
            Data_file = np.zeros((np.shape(index_05Z)[0], 150, 388))
            
            # loop through each day
            for dayx in range(len(index_05Z)):
                timex = arange(index_05Z[dayx], (index_05Z[dayx]+24))
                
                if (dayx == (len(index_05Z)-1)) & (t_mont[0] != '12'): 
                # pull following  month of data
                    first_range = arange(index_05Z[dayx], t_len)
                    first_halfv = getvar(ncfile,interp_method, timeidx=first_range, meta=False)
                    first_halfd = getvar(ncfile, fourd_var, timeidx=first_range, meta=False)
                    
                    second_halfv, second_halfd = secondFile(file_path, d03Files[filex+1], starting_hour, yr, interp_method, fourd_var)
                                        
                    vert       = mean(concatenate((first_halfv, second_halfv), axis = 0), axis=0)
                    var_windsx = mean(concatenate((first_halfd, second_halfd), axis = 0), axis=0)
                    var_winds = getvar(ncfile,fourd_var, timeidx=0, meta=True)
                    var_winds.data = var_windsx
                    
                    del first_range, first_halfv, first_halfd, second_halfv, second_halfd
                    
                    cs = vertcross(var_winds, vert, wrfin=ncfile, start_point=start_point, end_point=end_point, latlon=True, autolevels=150)
                    
                    del vert, var_windsx, var_winds
                    
                    Z_file, Data_file = processCrossSection(cs, dayx, Z_file, Data_file)
                    
                    del cs
                    
                elif (dayx == (len(index_05Z)-1)) & (t_mont[0] == '12') & (t_year[0] != str(end_yr)): 
                # Pull following year and month of data
                    first_range = arange(index_05Z[dayx], t_len)
                    first_halfv = getvar(ncfile,interp_method, timeidx=first_range, meta=False)
                    first_halfd = getvar(ncfile, fourd_var, timeidx=first_range, meta=False)
                    
                    file_path_next = filelocation + str(yr+1) + "/output/"
    
                    files_next = os.listdir(file_path_next)
                    d03Files_next = sorted([i for i in files_next if "wrfout_d03_" in i])
                    d03Files_next = [i for i in d03Files_next if "." not in i] # remove any bad files

                    second_halfv, second_halfd = secondFile(file_path_next, d03Files_next[0], starting_hour, yr+1, interp_method, fourd_var)
                    
                    vert       = mean(concatenate((first_halfv, second_halfv), axis = 0), axis=0)
                    var_windsx = mean(concatenate((first_halfd, second_halfd), axis = 0), axis=0)
                    var_winds = getvar(ncfile,fourd_var, timeidx=0, meta=True)
                    var_winds.data = var_windsx
                    
                    del first_range, first_halfv, first_halfd, file_path_next, d03Files_next, second_halfv, second_halfd
                    
                    cs = vertcross(var_winds, vert, wrfin=ncfile, start_point=start_point, end_point=end_point, latlon=True, autolevels=150)
                    
                    del vert, var_windsx, var_winds
                    
                    Z_file, Data_file = processCrossSection(cs, dayx, Z_file, Data_file)
                    
                    LATLON = [pair.latlon_str(fmt="{:.4f},{:.4f}") for pair in to_np(cs.coords["xy_loc"])]
                    LATLON = np.array([string.split(",") for string in LATLON])
                    LAT = LATLON[:,0]
                    LON = LATLON[:,1]
                    
                    del cs, LATLON
                    
                elif (dayx == (len(index_05Z)-1)) & (t_mont[0] == '12') & (t_year[0] == str(end_yr)):
                # pull the last few hour of data
                    first_range = arange(index_05Z[dayx], t_len)
                    
                    vert       = mean(getvar(ncfile,interp_method, timeidx=first_range, meta=False), axis=0)
                    var_windsx = mean(getvar(ncfile, fourd_var, timeidx=first_range, meta=False), axis=0)
                    var_winds = getvar(ncfile,fourd_var, timeidx=0, meta=True)
                    var_winds.data = var_windsx
                    
                    cs = vertcross(var_winds, vert, wrfin=ncfile, start_point=start_point, end_point=end_point, latlon=True, autolevels=150)
                    
                    del vert, var_windsx, var_winds, first_range
                    
                    Z_file, Data_file = processCrossSection(cs, dayx, Z_file, Data_file)
                    
                    LATLON = [pair.latlon_str(fmt="{:.4f},{:.4f}") for pair in to_np(cs.coords["xy_loc"])]
                    LATLON = np.array([string.split(",") for string in LATLON])
                    LAT = LATLON[:,0]
                    LON = LATLON[:,1]
                    
                    del cs, LATLON
                    
                else:
                    #open the vertical coordinate file   
                    var_winds = getvar(ncfile,fourd_var, timeidx=0, meta=True)
                    
                    vert       = mean(getvar(ncfile,interp_method, timeidx=timex, meta=False), axis=0)
                    var_windsx  = mean(getvar(ncfile, fourd_var, timeidx=timex, meta=False), axis=0)
                    var_winds.data = var_windsx
                    
                    cs = vertcross(var_winds, vert, wrfin=ncfile, start_point=start_point, end_point=end_point, latlon=True, autolevels=150)
                    
                    del vert, var_windsx, var_winds
                    
                    Z_file, Data_file = processCrossSection(cs, dayx, Z_file, Data_file)
                    
                    del cs

            ncfile.close()
            
            timecountend += len(t)
            
            Z[timecount:timecountend, :] = Z_file
            Data[timecount:timecountend, :, :] = Data_file
            
            timecount += len(t)
            
            del t, index_05Z, t_hour, t_year
            
        matfilename = 'wrf_' + str(fnum) + '_' + fourd_var + '_mer_cs_' + str(yr)
        
        LAT = LAT.astype(np.float)
        LON = LON.astype(np.float)
        
        LAT = np.ascontiguousarray(LAT)
        LON = np.ascontiguousarray(LON)
        Z = np.ascontiguousarray(Z)
        Data = np.ascontiguousarray(Data)

        data = {'LAT':LAT, 'LON':LON}
        
        # add vertical levels
        data.update({'Z':Z})
    
        # add data
        data.update({fourd_var[0]:Data[:,:,:]})
    
        #save data dictonary to .mat file
        sio.savemat(savelocation + matfilename + '.mat', data)
        
        del LAT, LON, Z, Data, data, Z_file, Data_file
    


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


def coordinates(latx, lonx):
    point = CoordPair(lat=latx, lon=lonx)
    
    return point


def secondFile(path, file, start_hour, year, interp_method, fourd_var):
    ncfile2 = nc.Dataset(path + file, 'r')
    t_hour2, t_year2, _, _, _ = pullYearHour(ncfile2)
    
    index_05Z2 = pullIndex05(start_hour, year, t_hour2, t_year2)
    
    second_range = arange(0, index_05Z2[0])
    second_vert = getvar(ncfile2,interp_method, timeidx=second_range, meta=False)
    second_var = getvar(ncfile2, fourd_var, timeidx=second_range, meta=False)
    
    ncfile2.close()

    return second_vert, second_var


def processCrossSection(cs, dayx, Z, Data):
    datax = cs.data
    Data[dayx,:,:] = datax

    Zx = to_np(cs.vertical)
    Z[dayx,:] = Zx
    
    return Z, Data



main()









