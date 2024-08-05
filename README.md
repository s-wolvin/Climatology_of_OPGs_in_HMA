# Climatology of OPGs in HMA
This repository stores the most useful programs created during my Masters Thesis Research at the University of Utah, funded by NASA's High Mountian Asia Team (HiMAT; Award 80NSSC20K1594). <br /> <br />

# Orographic Precipitation Gradient Downscaling Summary

Within High Mountain Asia (HMA), the annual melting of glaciers and snowpack provides vital freshwater to populations living downstream. Precipitation over HMA can directly affect the freshwater availability in this region by altering the mass balance of glaciers and snowpack. However, available reanalyses and downscaling simulations lack the resolution required to understand important glacier-scale variations in precipitation. This study aimed to determine the current characteristics of Orographic Precipitation Gradients (OPG) by curve-fitting daily precipitation as a function of elevation from a 15-year, 4-km grid spaced Weather Research and Forecasting (WRF) model simulation focused on the Himalayan, Karakoram, and Hindu-Kush mountain ranges. To facilitate precipitation curve-fitting, the WRF model grid points were separated into regions of similar orientation, referred to as facets. An Akaike Information Criterion-corrected value, an index to compare fit-model strengths, identified the need for a curvature term to account for a varying OPG with elevation. Regions with similar seasonal variability were found using $k$-means clustering of the monthly mean OPG coefficients. We evaluated the seasonal and intra-seasonal variability of the central Himalayan slope's OPG since most observational studies of OPG over HMA focus on this region. Initial testing of precipitation predictions using monthly coefficients showed promising results in downscaling daily WRF precipitation, with the daily mean absolute error at each grid point ranging from 0 to 12 mm. Results provide a physically-based context for machine learning algorithms being developed to predict OPG and downscale precipitation output from global climate models over HMA.


# Getting Started

 Many of these scripts require the following datasets/open-source code.

* Country Borders from Natural Earth: <br />
https://www.naturalearthdata.com/http//www.naturalearthdata.com/download/50m/cultural/ne_50m_admin_0_countries.zip <br />
* Export High-Resolution and Transparent Figures: <br />
https://www.mathworks.com/matlabcentral/fileexchange/23629-export_fig <br />
* M_MAP: A Mapping Package for MATLAB: <br />
https://www.eoas.ubc.ca/~rich/map.html


# Folder Structure
    .
    ├── correlation             #
    ├── curveFitting            # Evaluate, visualize, and quantify the relationship between precipitation and elevation.
    ├── facetingAlgorithm       # Group terrain into facets based on orientation. (Gibson et al. 1997)
    ├── kMeansClustering        # K-Means clustering of OPG coefficients
    ├── nclColormaps            # NCL Colormaps converted into MATLAB files.
    ├── opgVariability          # Evaluation of the seasonal variability on the Central Himalayan Slope
    ├── post-processingWRF      # Post-processing the WRF output files
    └── README.md                 
