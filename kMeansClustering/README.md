# kMeansClustering
### k_means_clustering.m
This MATLAB script pulls the OPG coefficients, z-scores each facet's OPGs based on their time series, and applies a k-means clustering algorithm. The k-means algorithm loops through 1 to the maximum set clusters and saves the results. 

### kmeans_coef_spread.m
This MATLAB script pulls the OPG coefficients, z-scores each facet's OPGs based on their time series, and applies a k-means clustering algorithm based on the number of clusters chosen and a list of starting random seeds. This produces a spread of k-means clustering outcomes, which are then matched up with one another to plot the distribution of OPG coefficient values based on the random seed. 

### p3_vs_avg_precip.m
This MATLAB script loads the k-means clustering outcome of k_means_clustering.m and plots the third coefficient (P3) of the defined cluster number next to the corresponding average precipitation to visualize the strong relationship between the two. 

### relabel_clusters.m
This MATLAB script reorders the clusters so that the clusters are labeled from largest to smallest. 
