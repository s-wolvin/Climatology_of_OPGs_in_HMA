# correlation
### finddatapoint.f90
This Fortran script, once compiled, it is a function used in **three_dim_corr_atmos_VS_opgCoef.py**. It loops through timeseries of every grid point and calculates its correlation coefficient to a timeseries of an OPG coefficient. For more information on using Fortran scripts as functions for Python scripts, go to: https://sites.engineering.ucsb.edu/~shell/che210d/f2py.pdf

### three_dim_corr_atmos_VS_opgCoef.py
This Python script pulls the OPG coefficients, an atmospheric variable at numerous pressure levels from the post-processed WRF output, and facet data. Then it calculates the correlation coefficient for every atmospheric grid point and single OPG coefficient timeseries'. It picks the location of highest correlation and calculates the linear regression slope and R^2 value. Then saves a .mat file containing the location of highest correlation, the correlation coefficient, linear regression slope, and R^2 value. 

