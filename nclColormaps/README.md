# NCL Colormaps
#### This folder contains many of the NCL colormaps found online, converted into MATLAB code. The entire colormap can be accessed by using the function name. <br />
*colormapName()*  <br /> <br />
#### The starting index, ending index, step size, whether the colormap is revesed, the total values, etc. can be called upon by the function. For example: <br />
*colormapName*('Start', 40) -- Include every color starting at index 40.

*colormapName*('End', 100, 'Reverse', true) -- Pull the first 100 colors, then reverse the order of the colors.

*colormapName*('Step', 2) -- Pull every-other color.

*colormapName*('Total', 10) -- Pull 10 total colors, evenly spaced within the list of colors. <br />

#### Total List of Arguments:
Reverse - Boolean:      Value to Indicate if Colormap will be Reversed <br />
Start   - Integer:      Start Index Value <br />
End     - Integer:      End Index Value <br />
Skip    - Int/Array:    Index Values to Skip <br />
Step    - Integer:      Step Index Value <br />
Total   - Integer:      Total Number of Colors Evenly Distributed <br />
Repeat  - Integer:      Number of times to repeat a color value <br />

## A Full List of NCL Colormaps Can Be Accessed At: <br />
https://www.ncl.ucar.edu/Document/Graphics/color_table_gallery.shtml
