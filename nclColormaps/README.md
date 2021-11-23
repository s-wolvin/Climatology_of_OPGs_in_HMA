# NCL Colormaps
#### This folder contains many of the NCL colormaps found online, converted into MATLAB code. The entire colormap can be accessed by using the function name. <br />
*colormapName()*  <br /> <br />
#### In addition, the starting index, ending index, step size, whether the colormap is revesed, and the total values can be called upon by the function. For example: <br />
*colormapName*('Start', 40) -- Include every color starting at index 40.

*colormapName*('End', 100, 'Reverse', true) -- Pull the first 100 colors, then reverse the order of the colors.

*colormapName*('Step', 2) -- Pull every-other color.

*colormapName*('Total', 10) -- Pull 10 total colors, evenly spaced within the list of colors. <br />

## A full list of NCL colormaps can be accessed at: <br />
https://www.ncl.ucar.edu/Document/Graphics/color_table_gallery.shtml
