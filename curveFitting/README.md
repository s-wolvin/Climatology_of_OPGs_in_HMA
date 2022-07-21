# curveFitting
### calc_himat_opg_quadratic_linear.m
This MATLAB script pulls the precipitation, elevation, and facet data. Fits both a linear and quadratic model to the precipitation-elevation gradinet on each facet. Then saves the fit-coefficients and statistics connected to each fit. These statistics include the p-values to the quadratic coefficients, the AICc and R^2 values for both the quadratic and linear fits. 

### plot_stats_quadratic_linear_OPG.m
This MATLAB scipts takes the statistical output files from **calc_himat_opg_quadratic_linear.m** to asscess the difference between the use of a linear or quadratic model for curve-fitting the OPG. It plots the percentage of precipitation events where the p-value is significant for the 1st coefficient of the quadratic fit, where the AICc_quadratic index is lower or equal to the AICc_linear index, and where the R^2 value is higher using a quadratic fit. Last, it plots the average percent difference between the R^2 values of the precipitation events. 
