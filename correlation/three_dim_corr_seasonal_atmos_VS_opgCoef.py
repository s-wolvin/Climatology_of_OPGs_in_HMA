#%% Savanna Wolvin
# Created: Jul 22nd, 2021
# Edited: Jul 22nd, 2022
    
# SUMMARY
# This script pulls an OPG coefficient and an atmospheric variable from the 
# post-processed WRF output and finds the grid point of highest correlation 
# for each season (i.e., DJF, MAM, JJA, SON).

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
# 'three_dim_lin_regress_season_' + var_name[poly_var_num] + '_' + wrf_nest + '_' + prdct_name
# .mat file which saves the grid point location of the highest correlation 
# between the OPG coefficient and the atmospheric variable (by season), along 
# with the linear regression slope, R^2 value, and correlation coefficient




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
poly_var_num = 1
var_name = ['p1','p2','p3']

# Atmospheric variable
prdct_name = 'wa'
prdct_vars = 'w'

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
    mat_file = sio.loadmat(opg_dir + opg_fname) # Load file
    
    tx = mat_file['t'] # Pull time
    idx2 = idx + len(tx)
    
    poly_var[:,idx:idx2,:] = mat_file['opg_polyfit2'] # Pull OPG coefficients
    
    idx = idx2

poly_var = poly_var[poly_var_num,:,:] # Pull only one OPG coefficient




#%% pull atmos var
print('Pull Atmospheric Variable...')

atmos_data = np.zeros((lay_size, np.shape(tref)[0], lat_size,lon_size))
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






#%% create array of indicies of each season
print("Formulate seasonal vectors...")

t_mont = tref.strftime('%m')
t_mont = t_mont.to_list()

djf = [idx for idx, elem in enumerate(t_mont) if (elem == '01' or elem == '02' or elem == '12')]
mam = [idx for idx, elem in enumerate(t_mont) if (elem == '03' or elem == '04' or elem == '05')]
jja = [idx for idx, elem in enumerate(t_mont) if (elem == '06' or elem == '07' or elem == '08')]
son = [idx for idx, elem in enumerate(t_mont) if (elem == '09' or elem == '10' or elem == '11')]




#%% calculate the correlations and the locations for each facet

### create size variable for fortran module 
lat_size = np.shape(atmos_data)[2]
lon_size = np.shape(atmos_data)[3]
z_size = np.shape(atmos_data)[0]
fi_size = np.shape(poly_var)[1]


### choose between new variables or to load up a previous variable
location = np.zeros((fi_size, 3, 4))
corrCoef = np.zeros((fi_size, 4))
linReg   = np.zeros((fi_size, 2, 4))
r_2      = np.zeros((fi_size, 4))


# mat_file = sio.loadmat(save_dir + 'three_dim_lin_regress_season_' + var_name[poly_var_num] + '_' + wrf_nest + '_' + prdct_name )
# location = mat_file['location']
# corrCoef = mat_file['corrCoef']
# linReg   = mat_file['linReg']
# r_2      = mat_file['r_2']


### loop through each facet
print("Begin looping through facets...")
for fi in np.arange(1,995+1):
# for fi in helpme:
    print(str(fi))
    
    ### if there is more than one datapoint on the facet, then continue
    if np.count_nonzero(facets == fi) > 2:        
        season = 0
        
        ### loop through each season
        for sx in [djf, mam, jja, son]:
    
            real_values = np.invert(np.isnan(poly_var[:,fi-1]))
            
            ### create array of precip & season days
            sx_days = np.zeros((len(tref)), dtype=bool) 
            sx_days[sx] = True
            timex = np.zeros((len(tref)), dtype=bool) 
                          
            for idx in range(len(real_values)-1): 
                # pull dates with real values for the OPG coefficients and 
                # dates within that season
                timex[idx] = sx_days[idx] & real_values[idx] 
                
            t_size = np.count_nonzero(timex)

            ### if there is more than one precipitation day in that season, 
            ### then continue
            if t_size > 1:
        
                ### calculate correlation and location
                corrx = finddatapoint.corr_coef_point(fi, \
                            atmos_data[:,timex[:],:,:], poly_var[timex,:], \
                                lon_size, lat_size, z_size, t_size, fi_size)
              
                ### remove infs
                corrx[abs(corrx) == float('inf')] = float('nan')
                
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
                location[fi-1, :, season] = idx + 1  
                    
                ### calculate linear regression
                alx = poly_var[timex, fi-1]

                px1 = atmos_data[idx[0], timex, idx[1], idx[2]]
                px1 = px1.reshape((-1,1))
                
                reg = LinearRegression().fit(px1,alx)

                r_2[fi-1, season] = reg.score(px1, alx)
                linReg[fi-1, 0, season] = reg.coef_
                linReg[fi-1, 1, season] = reg.intercept_
                    
            season += 1

    # save/update .mat array everytime a facet's correlation is calculated
    matfilename = 'three_dim_lin_regress_season_' + var_name[poly_var_num] + '_' + wrf_nest + '_' + prdct_name 
    data = {'location':location,'linReg':linReg,'r_2':r_2,'corrCoef':corrCoef}                 
    sio.savemat(save_dir + matfilename + '.mat', data)                    

    

