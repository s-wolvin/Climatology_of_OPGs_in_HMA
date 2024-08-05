"""
Savanna Wolvin
Created: May 30th, 2023
Edited: Jul 27th, 2023


##### SUMMARY #################################################################
For a single facet, plot the average atmosphere of a strong sublinear increase 
and superliner increase/sublinear decrease between DJF and JJA. Plots include 
horizontal plot of 500 hPa winds, 850 hPa heights and winds, and IVT, and a 
meridional cross-section plot of meridional moisture transport and transport 
anomalies.


"""
#%% Global Imports

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import cartopy.crs as ccrs
import cartopy.feature as cfeat
import scipy.io as sio
from datetime import datetime, timedelta
import scipy.stats as sps
import scipy as sp
import mat73 as mat
from tqdm import tqdm
# import xarray as xr
# import os
import sys
sys.path.insert(1, '/uufs/chpc.utah.edu/common/home/u1324060/nclcmappy/')
import nclcmaps as ncm


#%% Variable Presets

facet_dir   = "/uufs/chpc.utah.edu/common/home/strong-group4/savanna/wrf_facet_data/"
precip_dir  = "/uufs/chpc.utah.edu/common/home/strong-group4/savanna/precip_wrf/"
save_dir    = "/uufs/chpc.utah.edu/common/home/u1324060/himat_ms/pub_fig/"
poly_dir    = "/uufs/chpc.utah.edu/common/home/strong-group4/savanna/himat_opg_poly/"
wrf_d03_dir     = "/uufs/chpc.utah.edu/common/home/strong-group4/savanna/var_wrf/"
wrf_d02_dir = "/uufs/chpc.utah.edu/common/home/strong-group7/savanna/himat/jgr_wrf_cf_compress/d02/"

fnum = 939

years = [2001,2015]

press = [1000, 950, 900, 800, 700, 600, 500, 400, 300, 200]

# press_lev_djf = 3 # [900, 800, 700, 600, 500, 400, 300, 200]
# press_lev_jja = 3

prThrshld = 1

wrf_nest = "d02"

season = "DJF"


#%% Create Needed Domain Data

# Create Array of All Desired Years
data_years = np.arange(years[0], years[1]+1).astype('int')

# Create Array of All Desired Days
data_days = np.arange(datetime(years[0],1,1), datetime(years[1]+1,1,1), timedelta(days=1)).astype(datetime)

# Create Empty Variables
poly_var = np.zeros((np.shape(data_days)[0], 2))
pr_var  = np.zeros((np.shape(data_days)[0], 1))


#%% Load Facet Data

print('Load Precip and Coefficients...')

idx = 0

# Load Facet Location
mat_file = sio.loadmat(facet_dir + 'wrf_facets_25')
facets  = mat_file['facets_labeled']
use = (facets==fnum)

# Load Precipitation Data
for yx in tqdm(data_years):
    
    idx_ = np.shape(np.arange(datetime(yx,1,1), datetime(yx+1,1,1), timedelta(days=1)))[0] + idx
    
    # Load Precipitation
    mat_file = sio.loadmat(f"{precip_dir}daily_precip_{yx}")
    pr = mat_file['pr']
    pr_var[idx:idx_,0] = np.nanmean(pr[:,use], axis=1)

    # Load Polyfit Variables
    mat_file = mat.loadmat(f"{poly_dir}25km/daily_opg_{yx}_25km_polyfit2.mat")
    opg_polyfit2 = mat_file['opg_polyfit2']
    
    for dayx in range(np.shape(opg_polyfit2)[1]):
        poly_var[idx+dayx,0] = sps.mode(opg_polyfit2[0,dayx,use], axis=None, keepdims=False).mode
        poly_var[idx+dayx,1] = sps.mode(opg_polyfit2[1,dayx,use], axis=None, keepdims=False).mode

    idx = idx_

# Pull Percentile Values of P1 and P2

p1_q33 = np.nanpercentile(poly_var[:,0], 33)
p1_q66 = np.nanpercentile(poly_var[:,0], 66)

p2_q33 = np.nanpercentile(poly_var[:,1], 33)
p2_q66 = np.nanpercentile(poly_var[:,1], 66)


#%% Pull JJA and DJF

print('Define the seasons...')

jja_days = np.zeros((np.shape(data_days)[0], 1))
djf_days = np.zeros((np.shape(data_days)[0], 1))

for dayx in range(np.shape(data_days)[0]):
    mx = data_days[dayx].month
    if (mx == 1) or (mx == 2) or (mx == 12):
        djf_days[dayx,:] = 1
    elif (mx == 6) or (mx == 7) or (mx == 8):
        jja_days[dayx,:] = 1



#%% Load Atmospheric Data

print('Load Horizontal Atmos...')

idx = 0


if wrf_nest == "d03":
    geopt   = np.zeros((2, np.shape(data_days)[0], 402, 525))
    IVT     = np.zeros((np.shape(data_days)[0], 402, 525))
    uwnds   = np.zeros((2, np.shape(data_days)[0], 402, 525))
    vwnds   = np.zeros((2, np.shape(data_days)[0], 402, 525))
    
    for yx in tqdm(data_years):
        
        idx_ = np.shape(np.arange(datetime(yx,1,1), datetime(yx+1,1,1), timedelta(days=1)))[0] + idx
        
        # Load Geopotential Data
        mat_file = sio.loadmat(wrf_d03_dir + 'wrfout_d03_geopt_850_700_500_300mb_' + str(yx))
        geopt[0,idx:idx_,:,:] = mat_file['g'][0,:,:,:] # 850 hPa
        geopt[1,idx:idx_,:,:] = mat_file['g'][2,:,:,:] # 500 hPa
        
        # Load IVT
        mat_file = sio.loadmat(wrf_d03_dir + 'wrfout_d03_IVT_' + str(yx) + "_mk2")
        IVT[idx:idx_,:,:] = mat_file['IVT']
        
        
        # Load U-Winds Data
        mat_file = sio.loadmat(wrf_d03_dir + 'wrfout_d03_ua_1000_925_850mb_' + str(yx))
        uwnds[0,idx:idx_,:,:] = mat_file['u'][2,:,:,:] # 850 hPa
        
        mat_file = sio.loadmat(wrf_d03_dir + 'wrfout_d03_ua_700_600_500mb_' + str(yx))
        uwnds[1,idx:idx_,:,:] = mat_file['u'][2,:,:,:] # 500 hPa
        
        # Load V-Winds Data
        mat_file = sio.loadmat(wrf_d03_dir + 'wrfout_d03_va_1000_925_850mb_' + str(yx))
        vwnds[0,idx:idx_,:,:] = mat_file['v'][2,:,:,:] # 850 hPa
        
        mat_file = sio.loadmat(wrf_d03_dir + 'wrfout_d03_va_700_600_500mb_' + str(yx))
        vwnds[1,idx:idx_,:,:] = mat_file['v'][2,:,:,:] # 500 hPa
    
        
        idx = idx_
        
elif wrf_nest == "d02":
        
    geopt   = np.zeros((2, np.shape(data_days)[0], 219, 270))
    IVT     = np.zeros((np.shape(data_days)[0], 219, 270))
    uwnds   = np.zeros((2, np.shape(data_days)[0], 219, 270))
    vwnds   = np.zeros((2, np.shape(data_days)[0], 219, 270))
    
    for yx in tqdm(data_years):
    
        idx_ = np.shape(np.arange(datetime(yx,1,1), datetime(yx+1,1,1), timedelta(days=1)))[0] + idx
        
        # Load Geopotential Data
        mat_file = sio.loadmat(wrf_d03_dir + 'wrfout_d02_geopt_850_500mb_' + str(yx))
        geopt[0,idx:idx_,:,:] = mat_file['g'][0,:,:,:] # 850 hPa
        geopt[1,idx:idx_,:,:] = mat_file['g'][1,:,:,:] # 500 hPa
        
        # Load IVT
        mat_file = sio.loadmat(wrf_d03_dir + 'wrfout_d02_IVT_' + str(yx) + "_mk2")
        IVT[idx:idx_,:,:] = mat_file['IVT']
        
        
        # Load U-Winds Data
        mat_file = sio.loadmat(wrf_d03_dir + 'wrfout_d02_ua_1000_925_850mb_' + str(yx))
        uwnds[0,idx:idx_,:,:] = mat_file['u'][2,:,:,:] # 850 hPa
        
        mat_file = sio.loadmat(wrf_d03_dir + 'wrfout_d02_ua_700_600_500mb_' + str(yx))
        uwnds[1,idx:idx_,:,:] = mat_file['u'][2,:,:,:] # 500 hPa
        
        # Load V-Winds Data
        mat_file = sio.loadmat(wrf_d03_dir + 'wrfout_d02_va_1000_925_850mb_' + str(yx))
        vwnds[0,idx:idx_,:,:] = mat_file['v'][2,:,:,:] # 850 hPa
        
        mat_file = sio.loadmat(wrf_d03_dir + 'wrfout_d02_va_700_600_500mb_' + str(yx))
        vwnds[1,idx:idx_,:,:] = mat_file['v'][2,:,:,:] # 500 hPa
    
        
        idx = idx_

    
lat = mat_file['LAT']
lon = mat_file['LON']
    
mat_file = sio.loadmat(wrf_d03_dir + 'wrfout_d03_va_700_600_500mb_' + str(yx))
lat_fi = mat_file['LAT']
lon_fi = mat_file['LON']


del mat_file


#%% Remove extreme values

print('Remove outliers...')

geopt[geopt > 1000000] = np.nan
uwnds[uwnds > 1000000] = np.nan
vwnds[vwnds > 1000000] = np.nan

geopt[geopt < -900] = np.nan
uwnds[uwnds < -900] = np.nan
vwnds[vwnds < -900] = np.nan

geopt = geopt / 9.81
geopt = geopt / 10
# IWV   = IWV / 9.81
IVT   = IVT / 9.81


#%% Load in vertical data

print('Load vertical data...')

atmos_dir = "/uufs/chpc.utah.edu/common/home/strong-group4/savanna/var_wrf/"
wrf_nest = "wrfout_d02"
corr_atmos = "QVAPOR"
atmos_press = "900_800_700_600_500_400_300_200mb"
atmos_press2 = "1000_950_900mb"

Q           = np.zeros((10, np.shape(data_days)[0], 219))
hor_vwnds   = np.zeros((10, np.shape(data_days)[0], 219))
terrain_d02     = np.zeros((np.shape(data_days)[0], 219))
terrain_d03     = np.zeros((np.shape(data_days)[0], 402))

use = facets==939

difference_array = np.absolute(lon - np.mean(lon_fi[facets==939]))
index_d02 = np.argwhere(difference_array == np.min(difference_array[:]))
index_d02[0,1] = 185

difference_array = np.absolute(lon_fi - np.mean(lon_fi[facets==939]))
index_d03 = np.argwhere(difference_array == np.min(difference_array[:]))

idx = 0

for yx in tqdm(data_years):
    
    idx_ = np.shape(np.arange(datetime(yx,1,1), datetime(yx+1,1,1), timedelta(days=1)))[0] + idx

    # Load moisture
    mat_file = sio.loadmat(atmos_dir + wrf_nest + "_" + corr_atmos + "_" + atmos_press + "_" + str(yx))
    Q[2:10,idx:idx_,:] = mat_file['Q'][:,:,:,index_d02[0,1]]  * 1000 # 900_800_700_600_500_400_300_200mb, kg/kg to g/kg

    mat_file = sio.loadmat(atmos_dir + wrf_nest + "_" + corr_atmos + "_" + atmos_press2 + "_" + str(yx))
    Q[0:2,idx:idx_,:] = mat_file['Q'][0:2,:,:,index_d02[0,1]]  * 1000 # 1000_950_900mb, kg/kg to g/kg

    # load v-winds
    mat_file = sio.loadmat(atmos_dir + wrf_nest + "_" + "va" + "_" + atmos_press + "_" + str(yx))
    hor_vwnds[2:10,idx:idx_,:] = mat_file['v'][:,:,:,index_d02[0,1]] # 900_800_700_600_500_400_300_200mb

    mat_file = sio.loadmat(atmos_dir + wrf_nest + "_" + "va" + "_" + atmos_press2 + "_" + str(yx))
    hor_vwnds[0:2,idx:idx_,:] = mat_file['v'][0:2,:,:,index_d02[0,1]] # 1000_950_900mb
    
    # Load elevation
    mat_file = mat.loadmat(atmos_dir + wrf_nest + "_PSFC_" + str(yx) + "_24h_05Z.mat")
    terrain_d02[idx:idx_,:] = mat_file['wrf_var'][:,:,index_d02[0,1]]
    
    # mat_file = mat.loadmat(atmos_dir + "wrfout_d03_PSFC_" + str(yx) + "_24h_05Z.mat")
    # terrain_d03[idx:idx_,:] = mat_file['wrf_var'][:,:,index_d03[0,1]]
    
    idx = idx_

del mat_file

print('Remove outliers...')
Q[Q > 6000000]                  = np.nan
hor_vwnds[hor_vwnds > 6000000]  = np.nan
transport = Q * hor_vwnds
terrain_d02 = terrain_d02 / 100
# terrain_d03 = terrain_d03 / 100

 #%% Define facet latitude extent

fi_lat_d03 = lat_fi[:, index_d03[0,1]]
fi_lat_d03 = fi_lat_d03[use[:, index_d03[0,1]]]

fi_idx = np.logical_and(lat[:, index_d02[0,1]] > np.min(fi_lat_d03), 
                        lat[:, index_d02[0,1]] < np.max(fi_lat_d03))



#%%% Pull Average Atmos Values and smooth with gaussian filter

print("Formulate average atmos and smooth...")

if wrf_nest == "wrfout_d03":
    step = 80
    smooth = 2
elif wrf_nest == "wrfout_d02":
    step = 31
    smooth = 1

def filter_hgt(array):

    array = sp.ndimage.gaussian_filter(array, sigma=np.nanstd(array)/smooth)
    
    return array


if season == "JJA":
    # JJA Concave Up Pattern
    idx = (jja_days[:,0] == 1) & (poly_var[:,0] > p1_q66) & (poly_var[:,1] < p2_q33)
    
    jja_ccu_gph850  = filter_hgt(np.nanmean(geopt[0, idx, :, :], axis = 0))
    jja_ccu_gph500  = filter_hgt(np.nanmean(geopt[1, idx, :, :], axis = 0))
    jja_ccu_IVT     = np.nanmean(IVT[idx, :, :], axis = 0)
    jja_ccu_uwnd    = np.nanmean(uwnds[0, idx, :, :], axis = 0)[::step, ::step]
    jja_ccu_vwnd    = np.nanmean(vwnds[0, idx, :, :], axis = 0)[::step, ::step]
    
    jja_ccu_transport   = np.nanmean(transport[:, idx, :], axis = 1)
    jja_ccu_hor_vwnd    = np.nanmean(hor_vwnds[:, idx, :], axis = 1)
    
    # JJA Strong Concave Down Pattern
    idx = (jja_days[:,0] == 1) & (poly_var[:,0] < p1_q33) & (poly_var[:,1] > p2_q66)
    
    jja_ccd_gph850  = filter_hgt(np.nanmean(geopt[0, idx, :, :], axis = 0))
    jja_ccd_gph500  = filter_hgt(np.nanmean(geopt[1, idx, :, :], axis = 0))
    jja_ccd_IVT     = np.nanmean(IVT[idx, :, :], axis = 0)
    jja_ccd_uwnd    = np.nanmean(uwnds[0, idx, :, :], axis = 0)[::step, ::step]
    jja_ccd_vwnd    = np.nanmean(vwnds[0, idx, :, :], axis = 0)[::step, ::step]
    
    jja_ccd_transport   = np.nanmean(transport[:, idx, :], axis = 1)
    jja_ccd_hor_vwnd    = np.nanmean(hor_vwnds[:, idx, :], axis = 1)
    
    # Seasonal Mean
    jja_terrain         = np.nanmean(terrain_d02[jja_days[:,0] == 1, :], axis = 0)
    jja_mean_transport  = np.nanmean(transport[:, jja_days[:,0] == 1, :], axis = 1)

if season == "DJF":
    # DJF Concave Up Pattern
    idx = (djf_days[:,0] == 1) & (poly_var[:,0] > p1_q66) & (poly_var[:,1] < p2_q33)
    
    djf_ccu_gph850  = filter_hgt(np.nanmean(geopt[0, idx, :, :], axis = 0))
    djf_ccu_gph500  = filter_hgt(np.nanmean(geopt[1, idx, :, :], axis = 0))
    djf_ccu_IVT     = np.nanmean(IVT[idx, :, :], axis = 0)
    djf_ccu_uwnd    = np.nanmean(uwnds[1, idx, :, :], axis = 0)[::step, ::step]
    djf_ccu_vwnd    = np.nanmean(vwnds[1, idx, :, :], axis = 0)[::step, ::step]
    
    djf_ccu_transport   = np.nanmean(transport[:, idx, :], axis = 1)
    djf_ccu_hor_vwnd    = np.nanmean(hor_vwnds[:, idx, :], axis = 1)
    
    # DJF Strong Concave Down Pattern
    idx = (djf_days[:,0] == 1) & (poly_var[:,0] < p1_q33) & (poly_var[:,1] > p2_q66)
    
    djf_ccd_gph850  = filter_hgt(np.nanmean(geopt[0, idx, :, :], axis = 0))
    djf_ccd_gph500  = filter_hgt(np.nanmean(geopt[1, idx, :, :], axis = 0))
    djf_ccd_IVT     = np.nanmean(IVT[idx, :, :], axis = 0)
    djf_ccd_uwnd    = np.nanmean(uwnds[1, idx, :, :], axis = 0)[::step, ::step]
    djf_ccd_vwnd    = np.nanmean(vwnds[1, idx, :, :], axis = 0)[::step, ::step]
    
    djf_ccd_transport   = np.nanmean(transport[:, idx, :], axis = 1)
    djf_ccd_hor_vwnd    = np.nanmean(hor_vwnds[:, idx, :], axis = 1)
    
    # Seasonal Mean
    djf_terrain         = np.nanmean(terrain_d02[djf_days[:,0] == 1, :], axis = 0)
    djf_mean_transport  = np.nanmean(transport[:, djf_days[:,0] == 1, :], axis = 1)


# terrain values
hma_lat = np.append(lat[0,index_d02[0,1]], lat[:,index_d02[0,1]])
hma_lat = np.append(hma_lat, lat[-1,index_d02[0,1]])
if season == "DJF":
    hma_ter = np.append(1100, djf_terrain)
    hma_ter = np.append(hma_ter, 1100)
elif season == "JJA":
    hma_ter = np.append(1100, jja_terrain)
    hma_ter = np.append(hma_ter, 1100)


# Facet 939 values
fi_lat = lat[:,index_d02[0,1]]
fi_lat = fi_lat[fi_idx]
fi_lat = np.append(fi_lat[0], fi_lat)
fi_lat = np.append(fi_lat, fi_lat[-1])
if season == "DJF":
    fi_ter = djf_terrain[fi_idx]
    fi_ter = np.append(1000, fi_ter)
    fi_ter = np.append(fi_ter, 1000)
elif season == "JJA":
    fi_ter = jja_terrain[fi_idx]
    fi_ter = np.append(1000, fi_ter)
    fi_ter = np.append(fi_ter, 1000)







#%% Plot Presets

datacrs = ccrs.PlateCarree()
projex = ccrs.LambertConformal(central_longitude=np.mean(lon))
y_axis = 2
x_axis = 2
# extent

if wrf_nest == "wrfout_d03":
    extent = [np.min(lon), np.max(lon), np.min(lat)-1, np.max(lat)]
elif wrf_nest == "wrfout_d02":
    extent = [np.min(lon)+2, np.max(lon)-2, np.min(lat)-2, np.max(lat)]


IVT_levs = np.arange(0, 700, 100)

slp_levs = np.arange(992, 1022, 4)

# levs_500 = np.arange(5460, 6120, 60)
# levs_850 = np.arange(1380, 1600, 10)

levs_500 = np.arange(546, 612, 12)
levs_850 = np.arange(138, 160, 2)

landcolor = ''
bordercolor = [125/255, 125/255, 125/255]

fntsz = 10
lblsz = 9
ln_wth = 1.5


# Set Colors

borders = [0.5, 0.5, 0.5]

c_500 = 'black'
c_850 = 'purple'
c_wind = 'black'

IVT_cmap = ncm.cmapRange('WhiteYellowOrangeRed', start=0, finish=190)
tra_cmap = ncm.cmap('MPL_BrBG')

f_939 = 'blue'


# quiver scale
s_500 = 300
s_850 = 150


# meshgrid for cross-sections
XX, YY = np.meshgrid(lat[:,index_d02[0,1]], press)

# vapor transport levels
t_levs = np.arange(-35, 40, 5)


# manual locations
m_loc_before = [(69, 25), (69, 20), (88, 35)]
m_loc_1 = []
for lonx, latx in m_loc_before:
    m_loc_1.append(ccrs.LambertConformal(central_longitude=np.mean(lon)).transform_point(lonx, latx, ccrs.PlateCarree()))

m_loc_before = [(90, 24), (77, 24), (71, 28)]
m_loc_2 = []
for lonx, latx in m_loc_before:
    m_loc_2.append(ccrs.LambertConformal(central_longitude=np.mean(lon)).transform_point(lonx, latx, ccrs.PlateCarree()))
    
m_loc_before = [(68, 25), (69, 29), (82, 36)]
m_loc_3 = []
for lonx, latx in m_loc_before:
    m_loc_3.append(ccrs.LambertConformal(central_longitude=np.mean(lon)).transform_point(lonx, latx, ccrs.PlateCarree()))
    
m_loc_before = [(79, 24), (70, 31), (64, 38)]
m_loc_4 = []
for lonx, latx in m_loc_before:
    m_loc_4.append(ccrs.LambertConformal(central_longitude=np.mean(lon)).transform_point(lonx, latx, ccrs.PlateCarree()))
    
m_loc_before = [(63, 20), (79, 26), (64, 39)]
m_loc_5 = []
for lonx, latx in m_loc_before:
    m_loc_5.append(ccrs.LambertConformal(central_longitude=np.mean(lon)).transform_point(lonx, latx, ccrs.PlateCarree()))

m_loc_before = [(63, 18), (73, 26), (88, 38)]
m_loc_6 = []
for lonx, latx in m_loc_before:
    m_loc_6.append(ccrs.LambertConformal(central_longitude=np.mean(lon)).transform_point(lonx, latx, ccrs.PlateCarree()))
    
    

#%% Plot map and cross section DJF

if season == 'DJF':
    fig, ax = plt.subplots(nrows=y_axis, ncols=x_axis, figsize=(6.5, 4.75),
                           subplot_kw={'projection': projex})
    
    ##### DJF Strong Concave-Down Pattern - Map
    
    # Add outline of facet
    ax[0,0].contour(lon_fi, lat_fi, facets==939, colors=f_939, levels=[0.75, 10], 
                    linewidths=ln_wth, transform=datacrs, zorder=5)
    
    geo_con = ax[0,0].contour(lon, lat, djf_ccd_gph850, colors=c_850,
                     levels=levs_850, transform=datacrs, linewidths=ln_wth, zorder=7)
    
    ax[0,0].clabel(geo_con, fmt="%.0f", manual=m_loc_2, 
                          fontsize = lblsz, inline = True, 
                          use_clabeltext = True, inline_spacing = 25)
    
    ax[0,0].contourf(lon, lat, djf_ccd_IVT, cmap=IVT_cmap,
                    levels=IVT_levs, extend='max', transform=datacrs, zorder=2)
    
    qu = ax[0,0].quiver(lon[::step, ::step], lat[::step, ::step], djf_ccd_uwnd, djf_ccd_vwnd, 
                   scale=s_500, pivot='mid', transform=datacrs, color=c_500, zorder=8)
    
    ax[0,0].quiverkey(qu, 0.05, -0.02, 20, '20 $\\frac{m}{s}$', 
                      labelpos='S', coordinates='figure') 
    # fig.text(0.0, 0.86, "500 hPa:")
    
    ax[0,0].plot(lon[:,index_d02[0,1]], lat[:,index_d02[0,1]], '--', color=[0.2, 0.2, 0.2], transform=datacrs)
    
    ax[0,0].add_feature(cfeat.LAND, facecolor='lightgray', zorder=1)
    ax[0,0].add_feature(cfeat.COASTLINE.with_scale(
        '110m'), edgecolor=bordercolor, zorder=3)
    ax[0,0].add_feature(cfeat.BORDERS.with_scale(
        '50m'), edgecolor=bordercolor, zorder=4)
    ax[0,0].set_extent(extent)
    gl = ax[0,0].gridlines(crs=ccrs.PlateCarree(), \
                      xlocs=[60,75,90], \
                      ylocs=[20,30,40], draw_labels=True, x_inline=False, y_inline=False)
    gl.top_labels = False
    gl.bottom_labels = False
    gl.right_labels = False
    gl.ylabel_style = {'rotation': 0}
    
    ax[0,0].set_title("a) DJF Strong Sublinear Increase", 
                      fontsize=fntsz)
    
    
    ###### DJF Strong Concave-Down Pattern - Cross-section
    
    ax[0,1].remove()
    ax[0,1] = fig.add_subplot(y_axis, x_axis, 2, projection='rectilinear')
    
    ax[0,1].contourf(XX, YY, djf_ccd_transport, cmap=tra_cmap, 
                     levels=t_levs, extend='both')
    
    xy = np.stack((hma_lat, hma_ter), axis=1)
    ax[0,1].add_patch(patches.Polygon(xy, color='gray'))
    
    
    # ax[0,1].plot(fi_lat, fi_ter, 'b')
    
    xy = np.stack((fi_lat, fi_ter), axis=1)
    ax[0,1].add_patch(patches.Polygon(xy, facecolor='b'))
    
    anom = ax[0,1].contour(XX, YY, djf_ccd_transport-djf_mean_transport, 
                    np.arange(-50, 60, 10), colors='k')
    
    ax[0,1].clabel(anom, fmt="%.0f", fontsize = lblsz, inline = True, 
                          use_clabeltext = True, inline_spacing = 25)
    
    ax[0,1].set_ylim(ax[0,1].get_ylim()[::-1])
    ax[0,1].set_xticklabels([])
    
    plt.grid(True)
    
    ax[0,1].set_ylim([1000, 200])
    ax[0,1].set_ylabel("Pressure (hPa)")
    
    ax[0,1].set_title("b) DJF Strong Sublinear Increase", 
                      fontsize=fntsz)
    
    
    ##### DJF Concave Up Pattern - Map
    
    ax[1,0].contour(lon_fi, lat_fi, facets==939, colors=f_939, levels=[0.75, 10], 
                    linewidths=ln_wth, transform=datacrs, zorder=5)
    
    geo_con = ax[1,0].contour(lon, lat, djf_ccu_gph850, colors=c_850,
                      levels=levs_850, transform=datacrs, linewidths=ln_wth, zorder=7)
    ax[1,0].clabel(geo_con, fmt="%.0f", manual=m_loc_4, 
                          fontsize = lblsz, inline = True, 
                          use_clabeltext = True, inline_spacing = 25)
    
    # ax[0,1].clabel(geo_con, fmt="%.0f", inline=True, fontsize=lblsz, inline_spacing=25)
    
    ivt = ax[1,0].contourf(lon, lat, djf_ccu_IVT, cmap=IVT_cmap,
                    levels=IVT_levs, extend='max', transform=datacrs, zorder=2)
    
    ax[1,0].quiver(lon[::step, ::step], lat[::step, ::step], djf_ccu_uwnd, djf_ccu_vwnd, 
                    scale=s_500, pivot='mid', transform=datacrs, color=c_500, zorder=8)
    
    ax[1,0].plot(lon[:,index_d02[0,1]], lat[:,index_d02[0,1]], '--', color=[0.2, 0.2, 0.2], transform=datacrs)
    
    ax[1,0].add_feature(cfeat.LAND, facecolor='lightgray', zorder=1)
    ax[1,0].add_feature(cfeat.COASTLINE.with_scale(
        '110m'), edgecolor=bordercolor, zorder=3)
    ax[1,0].add_feature(cfeat.BORDERS.with_scale(
        '50m'), edgecolor=bordercolor, zorder=4)
    ax[1,0].set_extent(extent)
    gl = ax[1,0].gridlines(crs=ccrs.PlateCarree(), \
                      xlocs=[60,75,90], \
                      ylocs=[20,30,40], draw_labels=True, x_inline=False, y_inline=False)
    gl.top_labels = False
    # gl.bottom_labels = False
    gl.right_labels = False
    gl.xlabel_style = {'rotation': 35}
    gl.ylabel_style = {'rotation': 0}
    
    ax[1,0].set_title("c) DJF Superlinear Increase", fontsize=fntsz)
    
    cbar_ivt = fig.add_axes([0.12, -0.03, 0.35, 0.04])
    cb = plt.colorbar(ivt, cax=cbar_ivt, orientation='horizontal')
    cb.set_label('IVT (kg m$^{-1}$ s$^{-1}$)')
    
    ###### DJF Strong Concave-Up Pattern - Cross-section
    
    ax[1,1].remove()
    ax[1,1] = fig.add_subplot(y_axis, x_axis, 4, projection='rectilinear')
    
    mt = ax[1,1].contourf(XX, YY, djf_ccu_transport, cmap=tra_cmap, 
                          levels=t_levs, extend='both')
    
    plt.grid(True)
    
    xy = np.stack((hma_lat, hma_ter), axis=1)
    ax[1,1].add_patch(patches.Polygon(xy, color='gray'))
    
    xy = np.stack((fi_lat, fi_ter), axis=1)
    ax[1,1].add_patch(patches.Polygon(xy, facecolor='b'))
    
    anom = ax[1,1].contour(XX, YY, djf_ccu_transport-djf_mean_transport,
                    np.arange(-50, 60, 10), colors='k')
    
    ax[1,1].clabel(anom, fmt="%.0f", fontsize = lblsz, inline = True, 
                          use_clabeltext = True, inline_spacing = 25)
    
    ax[1,1].set_ylim(ax[0,1].get_ylim()[::-1])
    
    
    ax[1,1].set_ylim([1000, 200])
    ax[1,1].set_xticks([20,25,30,35,40])
    ax[1,1].set_xticklabels(['20$^o$N','25$^o$N','30$^o$N','35$^o$N','40$^o$N'])
    ax[1,1].set_ylabel("Pressure (hPa)")
    
    ax[1,1].set_title("d) DJF Superlinear Increase", fontsize=fntsz)
    
    cbar_mt = fig.add_axes([0.54, -0.03, 0.48, 0.04])
    cb = plt.colorbar(mt, cax=cbar_mt, orientation='horizontal')
    cb.set_label('Meridional Moisture Transport (g kg$^{-1}$ m s$^{-1}$)')
    
    fig.tight_layout(pad=0.6)
    
    save_dir = "/uufs/chpc.utah.edu/common/home/strong-group7/savanna/himat/new_pub_figs/"

    plt.savefig(f"{save_dir}map_{wrf_nest}_DJF_hrzntl_vrtcl_atmos.eps", dpi=600, 
                transparent=True, bbox_inches='tight')
    

elif season == "JJA":
    fig, ax = plt.subplots(nrows=y_axis, ncols=x_axis, figsize=(6.5, 4.75),
                           subplot_kw={'projection': projex})
    
    ##### JJA Strong Concave-Down Pattern - Map
    
    # Add outline of facet
    ax[0,0].contour(lon_fi, lat_fi, facets==939, colors=f_939, levels=[0.75, 10], 
                    linewidths=ln_wth, transform=datacrs, zorder=5)
    
    geo_con = ax[0,0].contour(lon, lat, jja_ccd_gph850, colors=c_850,
                     levels=levs_850, transform=datacrs, linewidths=ln_wth, zorder=7)
    
    ax[0,0].clabel(geo_con, fmt="%.0f", manual=m_loc_5, 
                          fontsize = lblsz, inline = True, 
                          use_clabeltext = True, inline_spacing = 25)
    
    ax[0,0].contourf(lon, lat, jja_ccd_IVT, cmap=IVT_cmap,
                    levels=IVT_levs, extend='max', transform=datacrs, zorder=2)
    
    qu = ax[0,0].quiver(lon[::step, ::step], lat[::step, ::step], jja_ccd_uwnd, jja_ccd_vwnd, 
                   scale=s_850, pivot='mid', transform=datacrs, color=c_850, zorder=8)
    
    ax[0,0].quiverkey(qu, 0.05, -0.02, 10, '10 $\\frac{m}{s}$', 
                      labelpos='S', coordinates='figure') 
    # fig.text(0.0, 0.86, "500 hPa:")
    
    ax[0,0].plot(lon[:,index_d02[0,1]], lat[:,index_d02[0,1]], '--', color=[0.2, 0.2, 0.2], transform=datacrs)
    
    ax[0,0].add_feature(cfeat.LAND, facecolor='lightgray', zorder=1)
    ax[0,0].add_feature(cfeat.COASTLINE.with_scale(
        '110m'), edgecolor=bordercolor, zorder=3)
    ax[0,0].add_feature(cfeat.BORDERS.with_scale(
        '50m'), edgecolor=bordercolor, zorder=4)
    ax[0,0].set_extent(extent)
    gl = ax[0,0].gridlines(crs=ccrs.PlateCarree(), \
                      xlocs=[60,75,90], \
                      ylocs=[20,30,40], draw_labels=True, x_inline=False, y_inline=False)
    gl.top_labels = False
    gl.bottom_labels = False
    gl.right_labels = False
    gl.ylabel_style = {'rotation': 0}
    
    ax[0,0].set_title("a) JJA Strong Sublinear Increase", 
                      fontsize=fntsz)
    
    
    ###### DJF Strong Concave-Down Pattern - Cross-section
    
    ax[0,1].remove()
    ax[0,1] = fig.add_subplot(y_axis, x_axis, 2, projection='rectilinear')
    
    ax[0,1].contourf(XX, YY, jja_ccd_transport, cmap=tra_cmap, 
                     levels=t_levs, extend='both')
    
    xy = np.stack((hma_lat, hma_ter), axis=1)
    ax[0,1].add_patch(patches.Polygon(xy, color='gray'))
    
    xy = np.stack((fi_lat, fi_ter), axis=1)
    ax[0,1].add_patch(patches.Polygon(xy, facecolor='b'))
    
    anom = ax[0,1].contour(XX, YY, jja_ccd_transport-jja_mean_transport, 
                    np.arange(-50, 60, 10), colors='k')
    
    ax[0,1].clabel(anom, fmt="%.0f", fontsize = lblsz, inline = True, 
                          use_clabeltext = True, inline_spacing = 25)
    
    ax[0,1].set_ylim(ax[0,1].get_ylim()[::-1])
    ax[0,1].set_xticklabels([])
    
    plt.grid(True)
    
    ax[0,1].set_ylim([1000, 200])
    ax[0,1].set_ylabel("Pressure (hPa)")
    
    ax[0,1].set_title("b) JJA Strong Sublinear Increase", 
                      fontsize=fntsz)
    
    
    ##### DJF Concave Up Pattern - Map    
    ax[1,0].contour(lon_fi, lat_fi, facets==939, colors=f_939, levels=[0.75, 10], 
                    linewidths=ln_wth, transform=datacrs, zorder=5)
    
    geo_con = ax[1,0].contour(lon, lat, jja_ccu_gph850, colors=c_850,
                      levels=levs_850, transform=datacrs, linewidths=ln_wth, zorder=7)
    ax[1,0].clabel(geo_con, fmt="%.0f", manual=m_loc_6, 
                          fontsize = lblsz, inline = True, 
                          use_clabeltext = True, inline_spacing = 25)
    
    # ax[0,1].clabel(geo_con, fmt="%.0f", inline=True, fontsize=lblsz, inline_spacing=25)
    
    ivt = ax[1,0].contourf(lon, lat, jja_ccu_IVT, cmap=IVT_cmap,
                    levels=IVT_levs, extend='max', transform=datacrs, zorder=2)
    
    ax[1,0].quiver(lon[::step, ::step], lat[::step, ::step], jja_ccu_uwnd, jja_ccu_vwnd, 
                    scale=s_850, pivot='mid', transform=datacrs, color=c_850, zorder=8)
    
    ax[1,0].plot(lon[:,index_d02[0,1]], lat[:,index_d02[0,1]], '--', color=[0.2, 0.2, 0.2], transform=datacrs)
    
    ax[1,0].add_feature(cfeat.LAND, facecolor='lightgray', zorder=1)
    ax[1,0].add_feature(cfeat.COASTLINE.with_scale(
        '110m'), edgecolor=bordercolor, zorder=3)
    ax[1,0].add_feature(cfeat.BORDERS.with_scale(
        '50m'), edgecolor=bordercolor, zorder=4)
    ax[1,0].set_extent(extent)
    gl = ax[1,0].gridlines(crs=ccrs.PlateCarree(), \
                      xlocs=[60,75,90], \
                      ylocs=[20,30,40], draw_labels=True, x_inline=False, y_inline=False)
    gl.top_labels = False
    # gl.bottom_labels = False
    gl.right_labels = False
    gl.xlabel_style = {'rotation': 35}
    gl.ylabel_style = {'rotation': 0}
    
    ax[1,0].set_title("c) JJA Sublinear Decrease", fontsize=fntsz)
    
    cbar_ivt = fig.add_axes([0.12, -0.03, 0.35, 0.04])
    cb = plt.colorbar(ivt, cax=cbar_ivt, orientation='horizontal')
    cb.set_label('IVT (kg m$^{-1}$ s$^{-1}$)')
    
    ###### DJF Strong Concave-Up Pattern - Cross-section
    
    ax[1,1].remove()
    ax[1,1] = fig.add_subplot(y_axis, x_axis, 4, projection='rectilinear')
    
    mt = ax[1,1].contourf(XX, YY, jja_ccu_transport, cmap=tra_cmap, 
                          levels=t_levs, extend='both')
    
    plt.grid(True)
    
    xy = np.stack((hma_lat, hma_ter), axis=1)
    ax[1,1].add_patch(patches.Polygon(xy, color='gray'))
    
    xy = np.stack((fi_lat, fi_ter), axis=1)
    ax[1,1].add_patch(patches.Polygon(xy, facecolor='b'))
    
    anom = ax[1,1].contour(XX, YY, jja_ccu_transport-jja_mean_transport,
                    np.arange(-50, 60, 10), colors='k')
    
    ax[1,1].clabel(anom, fmt="%.0f", fontsize = lblsz, inline = True, 
                          use_clabeltext = True, inline_spacing = 25)
    
    ax[1,1].set_ylim(ax[0,1].get_ylim()[::-1])
    
    
    ax[1,1].set_ylim([1000, 200])
    ax[1,1].set_xticks([20,25,30,35,40])
    ax[1,1].set_xticklabels(['20$^o$N','25$^o$N','30$^o$N','35$^o$N','40$^o$N'])
    ax[1,1].set_ylabel("Pressure (hPa)")
    
    ax[1,1].set_title("d) JJA Sublinear Decrease", fontsize=fntsz)
    
    cbar_mt = fig.add_axes([0.54, -0.03, 0.48, 0.04])
    cb = plt.colorbar(mt, cax=cbar_mt, orientation='horizontal')
    cb.set_label('Meridional Moisture Transport (g kg$^{-1}$ m s$^{-1}$)')
    
    fig.tight_layout(pad=0.6)
    # plt.subplots_adjust(right=0.95)
    
    save_dir = "/uufs/chpc.utah.edu/common/home/strong-group7/savanna/himat/new_pub_figs/"

    plt.savefig(f"{save_dir}map_{wrf_nest}_JJA_hrzntl_vrtcl_atmos.eps", dpi=600, 
                transparent=True, bbox_inches='tight')




#%% Plot Variables

fig, ax = plt.subplots(nrows=y_axis, ncols=x_axis, figsize=(6, 4.75),
                       subplot_kw={'projection': projex})

##### DJF Strong Concave-Down Pattern
geo_con = ax[0,0].contour(lon, lat, djf_ccd_gph500, colors=c_500,
                 levels=levs_500, transform=datacrs, linewidths=1, zorder=6)

ax[0,0].clabel(geo_con, fmt="%.0f", manual=m_loc_1, 
                     fontsize = lblsz, inline = True, 
                     use_clabeltext = True, inline_spacing = 25)

# Add outline of facet
ax[0,0].contour(lon_fi, lat_fi, facets==939, colors=f_939, levels=[0.75, 10], 
                linewidths=1, transform=datacrs, zorder=5)

geo_con = ax[0,0].contour(lon, lat, djf_ccd_gph850, colors=c_850,
                 levels=levs_850, transform=datacrs, linewidths=1, zorder=7)

ax[0,0].clabel(geo_con, fmt="%.0f", manual=m_loc_2, 
                     fontsize = lblsz, inline = True, 
                     use_clabeltext = True, inline_spacing = 25)

ax[0,0].contourf(lon, lat, djf_ccd_IVT, cmap=IVT_cmap,
                levels=IVT_levs, transform=datacrs, zorder=2)

qu = ax[0,0].quiver(lon[::step, ::step], lat[::step, ::step], djf_ccd_uwnd, djf_ccd_vwnd, 
               scale=s_500, pivot='mid', transform=datacrs, color=c_500, zorder=8)

ax[0,0].quiverkey(qu, 0.99, 0.83, 20, '20 $\\frac{m}{s}$', 
                  labelpos='W', coordinates='figure') 
fig.text(0.92, 0.86, "500 hPa:")

ax[0,0].add_feature(cfeat.LAND, facecolor='lightgray', zorder=1)
ax[0,0].add_feature(cfeat.COASTLINE.with_scale(
    '110m'), edgecolor=bordercolor, zorder=3)
ax[0,0].add_feature(cfeat.BORDERS.with_scale(
    '50m'), edgecolor=bordercolor, zorder=4)
ax[0,0].set_extent(extent)
gl = ax[0,0].gridlines(crs=ccrs.PlateCarree(), \
                  xlocs=[60,75,90], \
                  ylocs=[20,30,40], draw_labels=True, x_inline=False, y_inline=False)
gl.top_labels = False
gl.bottom_labels = False
gl.right_labels = False
gl.ylabel_style = {'rotation': 0}

ax[0,0].set_title("a) DJF Strong Sublinear Increase", 
                  fontsize=fntsz)


##### DJF Concave Up Pattern
geo_con = ax[0,1].contour(lon, lat, djf_ccu_gph500, colors=c_500, 
                 levels=levs_500, transform=datacrs, linewidths=1, zorder=6)
ax[0,1].clabel(geo_con, fmt="%.0f", manual=m_loc_3, 
                     fontsize = lblsz, inline = True, 
                     use_clabeltext = True, inline_spacing = 25)

# ax[0,1].clabel(geo_con, geo_con.levels[::2], fmt="%.0f", inline=True, fontsize=lblsz, inline_spacing=25)

ax[0,1].contour(lon_fi, lat_fi, facets==939, colors=f_939, levels=[0.75, 10], 
                linewidths=1, transform=datacrs, zorder=5)

geo_con = ax[0,1].contour(lon, lat, djf_ccu_gph850, colors=c_850,
                 levels=levs_850, transform=datacrs, linewidths=1, zorder=7)
ax[0,1].clabel(geo_con, fmt="%.0f", manual=m_loc_4, 
                     fontsize = lblsz, inline = True, 
                     use_clabeltext = True, inline_spacing = 25)

# ax[0,1].clabel(geo_con, fmt="%.0f", inline=True, fontsize=lblsz, inline_spacing=25)

ax[0,1].contourf(lon, lat, djf_ccu_IVT, cmap=IVT_cmap,
                levels=IVT_levs, transform=datacrs, zorder=2)

ax[0,1].quiver(lon[::step, ::step], lat[::step, ::step], djf_ccu_uwnd, djf_ccu_vwnd, 
               scale=s_500, pivot='mid', transform=datacrs, color=c_500, zorder=8)

ax[0,1].add_feature(cfeat.LAND, facecolor='lightgray', zorder=1)
ax[0,1].add_feature(cfeat.COASTLINE.with_scale(
    '110m'), edgecolor=bordercolor, zorder=3)
ax[0,1].add_feature(cfeat.BORDERS.with_scale(
    '50m'), edgecolor=bordercolor, zorder=4)
ax[0,1].set_extent(extent)
ax[0,1].gridlines(crs=ccrs.PlateCarree(), \
                  xlocs=[60,75,90], \
                  ylocs=[20,30,40])

ax[0,1].set_title("b) DJF Superlinear Increase", fontsize=fntsz)


##### JJA Strong Concave-Down Pattern
geo_con = ax[1,0].contour(lon, lat, jja_ccd_gph500, colors=c_500, 
                 levels=levs_500, transform=datacrs, linewidths=1, zorder=6)
ax[1,0].clabel(geo_con, fmt="%.0f", inline=True, fontsize=lblsz, inline_spacing=25)

ax[1,0].contour(lon_fi, lat_fi, facets==939, colors=f_939, levels=[0.75, 10], 
                linewidths=1, transform=datacrs, zorder=5)

geo_con = ax[1,0].contour(lon, lat, jja_ccd_gph850, colors=c_850,
                 levels=levs_850, transform=datacrs, linewidths=1, zorder=7)
ax[1,0].clabel(geo_con, fmt="%.0f", manual=m_loc_5, 
                     fontsize = lblsz, inline = True, 
                     use_clabeltext = True, inline_spacing = 25)

# ax[1,0].clabel(geo_con, fmt="%.0f", inline=True, fontsize=lblsz, inline_spacing=25)

ax[1,0].contourf(lon, lat, jja_ccd_IVT, cmap=IVT_cmap,
                levels=IVT_levs, transform=datacrs, zorder=2)

qx = ax[1,0].quiver(lon[::step, ::step], lat[::step, ::step], jja_ccd_uwnd, jja_ccd_vwnd, 
               scale=s_850, pivot='mid', transform=datacrs, color=c_850, zorder=8)


ax[1,0].quiverkey(qx, 0.99, 0.15, 10, r'10 $\frac{m}{s}$', 
                  labelpos='W', coordinates='figure')
fig.text(0.92, 0.18, "850 hPa:")

ax[1,0].add_feature(cfeat.LAND, facecolor='lightgray', zorder=1)
ax[1,0].add_feature(cfeat.COASTLINE.with_scale(
    '110m'), edgecolor=bordercolor, zorder=3)
ax[1,0].add_feature(cfeat.BORDERS.with_scale(
    '50m'), edgecolor=bordercolor, zorder=4)
ax[1,0].set_extent(extent)
gl = ax[1,0].gridlines(crs=ccrs.PlateCarree(), \
                  xlocs=[60,75,90], \
                  ylocs=[20,30,40], draw_labels=True, x_inline=False, y_inline=False)
gl.top_labels = False
# gl.bottom_labels = False
gl.right_labels = False
gl.xlabel_style = {'rotation': 45}
gl.ylabel_style = {'rotation': 0}

ax[1,0].set_title("c) JJA Strong Sublinear Increase", fontsize=fntsz)


##### JJA Concave Up Pattern
geo_con = ax[1,1].contour(lon, lat, jja_ccu_gph500, colors=c_500, 
                 levels=levs_500, transform=datacrs, linewidths=1, zorder=6)
ax[1,1].clabel(geo_con, fmt="%.0f", inline=True, fontsize=lblsz, inline_spacing=75)

ax[1,1].contour(lon_fi, lat_fi, facets==939, colors=f_939, levels=[0.75, 10], 
                linewidths=1, transform=datacrs, zorder=5)

geo_con = ax[1,1].contour(lon, lat, jja_ccu_gph850, colors=c_850,
                 levels=levs_850, transform=datacrs, linewidths=1, zorder=7)
ax[1,1].clabel(geo_con, fmt="%.0f", manual=m_loc_6, 
                     fontsize = lblsz, inline = True, 
                     use_clabeltext = True, inline_spacing = 25)

# ax[1,1].clabel(geo_con, fmt="%.0f", inline=True, fontsize=lblsz, inline_spacing=25)

geo_ax = ax[1,1].contourf(lon, lat, jja_ccu_IVT, cmap=IVT_cmap,
                levels=IVT_levs, transform=datacrs, zorder=2)

ax[1,1].quiver(lon[::step, ::step], lat[::step, ::step], jja_ccu_uwnd, jja_ccu_vwnd, 
               scale=s_850, pivot='mid', transform=datacrs, color=c_850, zorder=8)

ax[1,1].add_feature(cfeat.LAND, facecolor='lightgray', zorder=1)
ax[1,1].add_feature(cfeat.COASTLINE.with_scale(
    '110m'), edgecolor=bordercolor, zorder=3)
ax[1,1].add_feature(cfeat.BORDERS.with_scale(
    '50m'), edgecolor=bordercolor, zorder=4)
ax[1,1].set_extent(extent)
gl = ax[1,1].gridlines(crs=ccrs.PlateCarree(), \
                  xlocs=[60,75,90], \
                  ylocs=[20,30,40], draw_labels=True, x_inline=False, y_inline=False)
gl.top_labels = False
gl.bottom_labels = True
gl.right_labels = False
gl.left_labels = False
gl.xlabel_style = {'rotation': 45}

ax[1,1].set_title("d) JJA Sublinear Decrease", fontsize=fntsz)


cbar_atmos = fig.add_axes([0.92, 0.26, 0.03, 0.5]) # x0, y0, dx, dy
cb_at = fig.colorbar(geo_ax, cax=cbar_atmos,
                      pad=0.0, aspect=15, fraction=0.032)
cb_at.set_label('IVT (kg m-1 s-1)', size=10)
cb_at.ax.tick_params(labelsize=10)

plt.subplots_adjust(wspace=0.05, hspace=0.25)

save_dir = "/uufs/chpc.utah.edu/common/home/strong-group7/savanna/himat/new_pub_figs/"

plt.savefig(f"{save_dir}map_{wrf_nest}_IVT_850hpa_500hpa_hgts_wnds_mk2.eps", dpi=600, 
            transparent=True, bbox_inches='tight')


plt.show()





