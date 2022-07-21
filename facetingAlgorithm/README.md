# facetingAlgorithm
### gibson_tree.m
Function used by the **wrf_facets_MAIN.m** script. It is the algorithm formulated by Wayne Gibson, Christopher Daly and George Taylor from Oregon State University, Corvallis, Oregon. It is used by the Parameter-elevation Regressions on Independent Slopes Model (PRISM). The algoirthm outlined can be found at https://prism.oregonstate.edu/documents/pubs/1997cac_derivationGrids_daly.pdf

### orientation_is_adjacent.m
Function needed for **gibson_tree.m** algorithm. Checks if two orientations are adjacent to one another. 

### orientation_sum.m
Function needed for **gibson_tree.m** algorithm. Calculates sum of grid points at orientaions adjacent to one another. 

### wrf_facets_MAIN.m
MATLAB script which loads in terrain data from a WRF output file, plots the smoothed terrain, plots flat vs not flat terrain, runs each gridpoint through the faceting alorithm outlined in gibson_tree.m, plots the determined orientations of each gridpoint, and saves a MAT file containing arrays of the orientations and the new labeled facets. This algorithm is used to facet the terrain (i.e., turn terrain into areas grouped by similar facing orienation, like mountain faces). It allows for the grouping of grid points to be used to evaluate the relationship between elevation and any atmospheric variable.
