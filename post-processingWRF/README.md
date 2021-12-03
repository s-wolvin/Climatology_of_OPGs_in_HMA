# post-processingWRF
### pull_himat_wrfVars_24hrAvg.m
This MATLAB script pulls any variable with a standard grid from the HiMAT WRF Model, creates daily averages from any chosen start time, and then saves it as a .mat file. 

### pull_himat_wrfVars_1hr.m
This MATLAB script pulls any variable with a standard grid from the HiMAT WRF Model, isolates the daily observations of a specified UTC hour, and then saves it as a .mat file. 

### pull_plot_himat_minMax_tempFrequency_hourly.m
This MATLAB script pulls the temperature data from the HiMAT WRF Model, counts the occurences of minimum and maximum temperatures at each hour at each gridpoint, and then plots the frequency in a .png file. 

### pull_plot_himat_precipFrequency_hourly.m
This MATLAB script pulls the precipitation data from the HiMAT WRF Model, counts the occurences of precipiation at rates > 3 mm/hr, >1 mm/hr && < 3 mm/hr, and < 1 mm/hr at each gridpoint and calculates the sum of all precipiation that occurs at each hour. Then plots the frequency of differing precipitation rates at each hour and the percentage of total precipitation that occurs at each hour over the entire simulation as .png files. 
