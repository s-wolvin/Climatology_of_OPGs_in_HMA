#%% Savanna Wolvin
# Created: Jul 22nd, 2021
# Edited: Jul 22nd, 2022
    
# SUMMARY
# This script pulls an OPG coefficient and an atmospheric variable from the 
# post-processed WRF output and finds the grid point of highest correlation.

# To run this script, you must compile finddatapoint.f90, which is a fortran 
# script created to speed up the calculation of correlation, which can involve 
# millions of loops

# To compile finddatapoint.f90, run the following linux command:
    # f2py -c -m corr_coef_point finddatapoint.f90 --fcompiler=gnu95 --compiler=mingw32

# For more information on creating python modules using fortran:
    # https://sites.engineering.ucsb.edu/~shell/che210d/f2py.pdf

# INPUT
# radius -          Radius or radii used for faceting the terrain from 
#                   wrf_facets_MAIN.m
# start_yr -        Starting year of precipitation data to loop through
# end_yr -          Ending year of precipitation data to loop through
# poly_fit - 
# poly_var_num -    Coefficient number to use for correlation
# var_name -        Names of the coefficients
# prdct_name -      Name of atmospheric variable to use for correlation, i.e.,
#                   name on .mat file
# prdct_vars -      Atmospheric variable saved within the .mat file
# wrf_nest -        WRF nest to pull atmoshperic data
# lat_size -        Number of grid points latitudinally
# lon_size -        Number of grid points longitudinally
# lay_size -        Number of pressure levels
# starting_hour -   Starting hour of the 24-hour average atmospheric variables
# opg_dir -         Directory of the OPG coefficients
# prdct_dir -       Directory of the atmospheric variables
# facet_dir -       Directory of the facet data
# save_dir -        Directory to save to

# OUTPUT
# 'three_dim_lin_regress_' + var_name[poly_var_num] + '_' + wrf_nest + '_' + prdct_name
# .mat file which saves the grid point location of the highest correlation 
# between the OPG coefficient and the atmospheric variable, along with the 
# linear regression slope, R^2 value, and correlation coefficient




#%% Imports
import numpy as np
import pandas as pd
from datetime import timedelta, datetime
import scipy.io as sio
from sklearn.linear_model import LinearRegression
import os




#%% Variable presets

# Sampling radius of the faceting algorithm
radius = '25'

# years
start_yr = 2001
end_yr = 2015

# Which quadratic variables to evaluate
poly_fit = 'polyfit2'
poly_var_num = 0
var_name = ['p1','p2','p3']

# Atmospheric variable
prdct_name = 'va'
prdct_vars = 'v'

# WRF output information
wrf_nest = 'wrfout_d02'
lat_size = 219 # D03: 402, D02: 219
lon_size = 270 # D03: 525, D02: 270
lay_size = 8   # Number of pressure levels
starting_hour = '05'

# Check our location, pick correct file paths
cwd = os.getcwd()

if cwd[1:5] == 'uufs':
    print('Running Remotely')
    import finddatapoint
    
    opg_dir = '/uufs/chpc.utah.edu/common/home/strong-group4/savanna/himat_opg_poly/' + radius + 'km/'
    prdct_dir = '/uufs/chpc.utah.edu/common/home/strong-group4/savanna/var_wrf/'
    facet_dir = '/uufs/chpc.utah.edu/common/home/strong-group4/savanna/wrf_facet_data/'
    save_dir = '/uufs/chpc.utah.edu/common/home/strong-group4/savanna/himat_opg_poly/'
    
else:
    print('Running Locally')
    import sys
    sys.path.append('C:/Users/u1324060/Documents/thesis')  
    import finddatapoint
    
    opg_dir = 'Z:/himat_opg_poly/' + radius + 'km/'
    prdct_dir = 'Z:/var_wrf/'
    facet_dir = 'Z:/wrf_facet_data/'
    save_dir = 'Z:/himat_opg_poly/'




#%% load facets
print('Load Facets and Date Vector...')

# Load facets
facets_file = sio.loadmat(facet_dir + 'wrf_facets_' + radius)
facets = facets_file['facets_labeled']
num_of_facets = np.max(facets)
del facets_file

# create reference array of datetimes
years = np.arange(start_yr,end_yr+1)
tref = np.arange(datetime(start_yr, 1, 1, int(starting_hour)), \
                 datetime(end_yr+1, 1, 1, int(starting_hour)), \
                     timedelta(days=1), dtype='object')
    
tref = pd.to_datetime(tref)




#%% pull poly vars
print('Pull Quadratic Coefficients...')

poly_var = np.zeros((3, np.shape(tref)[0], num_of_facets))
idx = 0

for yr in years:
    print(str(yr))
    opg_fname = 'daily_opg_' + str(yr) + '_' + radius + 'km_' + poly_fit + '_faceted'
    mat_file = sio.loadmat(opg_dir + opg_fname) # Load File
    
    tx = mat_file['t'] # Pull Time
    idx2 = idx + len(tx)
    
    poly_var[:,idx:idx2,:] = mat_file['opg_polyfit2'] # Pull OPG coefficients
    
    idx = idx2

poly_var = poly_var[poly_var_num,:,:] # Pull only one OPG coefficient




#%% pull atmos var
print('Pull Atmosphic Variable...')

atmos_data = np.zeros((lay_size, np.shape(tref)[0], lat_size,lon_size), dtype=np.single)
idx = 0

if wrf_nest == 'wrfout_d03':
    for yr in years:
        print(str(yr))
        prdct_fname = wrf_nest + '_' + prdct_name + '_850_800_300_200mb_' + str(yr)
        mat_file = sio.loadmat(prdct_dir + prdct_fname) # Load File
        
        atm = mat_file[prdct_vars] # Pull Var
        idx2 = idx + np.shape(atm)[1] # Find Array Size
        del atm
        
        atmos_data[[0,1,6,7], idx:idx2,:,:] = mat_file[prdct_vars] # Save Var to Array
        
        prdct_fname = wrf_nest + '_' + prdct_name + '_700_600_500_400mb_' + str(yr)
        mat_file = sio.loadmat(prdct_dir + prdct_fname) # Load File
        
        atmos_data[[2,3,4,5], idx:idx2,:,:] = mat_file[prdct_vars] # Save Var to Array
        
        idx = idx2

    del mat_file
    
elif wrf_nest == 'wrfout_d02':
    for yr in years:
        print(str(yr))
        prdct_fname = wrf_nest + '_' + prdct_name + '_900_800_700_600_500_400_300_200mb_' + str(yr)
        mat_file = sio.loadmat(prdct_dir + prdct_fname) # Load File
        
        atm = mat_file[prdct_vars] # Pull Var
        idx2 = idx + np.shape(atm)[1] # Find Array Size
        del atm
        
        atmos_data[:, idx:idx2,:,:] = mat_file[prdct_vars] # Save Var to Array
        
        idx = idx2
    
    del mat_file
    
    
    

#%% calculate the correlations and the locations for each facet

### create size variable for fortran module 
lat_size = np.shape(atmos_data)[2]
lon_size = np.shape(atmos_data)[3]
z_size = np.shape(atmos_data)[0]
fi_size = np.shape(poly_var)[1]

### choose between new variables or to load up a previous variable
location = np.zeros((fi_size, 3))
corrCoef = np.zeros((fi_size, 1))
linReg   = np.zeros((fi_size, 2))
r_2      = np.zeros((fi_size, 1))
 
# mat_file = sio.loadmat(save_dir + 'three_dim_lin_regress_' + prdct_name )
# location = mat_file['location']
# corrCoef = mat_file['corrCoef']
# linReg   = mat_file['linReg']
# r_2      = mat_file['r_2']
# del mat_file


### loop through each facet
for fi in np.arange(1,995+1):
    print(str(fi))
    
    ### Check if facet has enough grid points
    if np.count_nonzero(facets == fi) > 2:
        real_values = np.invert(np.isnan(poly_var[:,fi-1]))
        t_size = np.count_nonzero(real_values)
        
        ### Check if facet recieved more than 1 day of precipitation
        if t_size > 1:
    
            ### Find grid point with highest correlation to quadratic OPG coefficient
            corrx = finddatapoint.corr_coef_point(fi, \
                            atmos_data[:,real_values,:,:], poly_var[real_values,:], \
                            lon_size, lat_size, z_size, t_size, fi_size)
                
            # Find location of highest correlation
            corrx = np.absolute(corrx)
            max_corrx = np.nanmax(corrx)
            idx = np.where(corrx == max_corrx)
            corrCoef[fi-1,0] = max_corrx
            
            ### Pull coordinates of highest correlation
            idx1 = idx[0]
            idx2 = idx[1]
            idx3 = idx[2]
            
            if len(idx1) == 1:
                idx = np.concatenate(idx)
            else:
                idx = np.array([idx1[0], idx2[0], idx3[0]])

            ### save location in MATLAB indexing
            location[fi-1, :] = idx + 1
            
            ### calculate linear regression
            alx = poly_var[real_values, fi-1]
                    
            px1 = atmos_data[idx[0], real_values, idx[1], idx[2]]
            px1 = px1.reshape((-1,1))
                    
            reg = LinearRegression().fit(px1,alx)
            
            r_2[fi-1, 0] = reg.score(px1, alx)
            linReg[fi-1, 0] = reg.coef_
            linReg[fi-1, 1] = reg.intercept_
                    

    matfilename = 'three_dim_lin_regress_' + var_name[poly_var_num] + '_' + wrf_nest + '_' + prdct_name
    data = {'location':location,'linReg':linReg,'r_2':r_2,'corrCoef':corrCoef}                 
    sio.savemat(save_dir + matfilename + '.mat', data)                    
    


