function [cmap] = amwg_blueyellowred(args)
% Edited: Nov. 11th, 2021
%
% Colormap Reference:
% https://www.ncl.ucar.edu/Document/Graphics/ColorTables/amwg_blueyellowred.shtml
%
% ncolors = 16
% cmap = amwg_blueyellowred()
% cmap = amwg_blueyellowred(argumentName, argumentValue)
% 
% INPUT ARGUMENTS:
%   Reverse - Boolean:      Value to Indicate if Colormap will be Reversed
%   Start   - Integer:      Start Index Value
%   End     - Integer:      End Index Value
%   Skip    - Int/Array:    Index Values to Skip
%   Step    - Integer:      Step Index Value
%   Total   - Integer:      Total Number of Colors Evenly Distributed
%   Repeat  - Integer:      Number of times to repeat a color value

arguments
    args.Start (1,1) {mustBeNumeric}
    args.End (1,1) {mustBeNumeric}
    args.Step (1,1) {mustBeNumeric}
    args.Total (1,1) {mustBeNumeric}
    args.Reverse
    args.Repeat
    args.Skip
end

cmap = [[130, 32, 240];...
    [0, 0, 150];...
    [0, 0, 205];...
    [65, 105, 225];...
    [30, 144, 255];...
    [0, 191, 255];...
    [160, 210, 255];...
    [210, 245, 255];...
    [255, 255, 200];...
    [255, 225, 50];...
    [255, 170, 0];...
    [255, 110, 0];...
    [255, 0, 0];...
    [200, 0, 0];...
    [160, 35, 35];...
    [255, 105, 180]];

if isfield(args, 'Reverse')
    if args.Reverse == true
        cmap = flipud(cmap);
    end
end

if isfield(args, 'Start')
    if isfield(args, 'End')
        if isfield(args, 'Skip')
            idx = args.Start:args.End;
            mem = ismember(idx, args.Skip);
            cmap = cmap(idx(~mem), :);
        else
            cmap = cmap(args.Start:args.End, :);
        end  
    elseif isfield(args, 'Skip')
        idx = args.Start:size(cmap,1);
        mem = ismember(idx, args.Skip);
        cmap = cmap(idx(~mem), :);
    else
        cmap = cmap(args.Start:end, :);
    end    
elseif isfield(args, 'End')
    if isfield(args, 'Skip')
        idx = 1:args.End;
        mem = ismember(idx, args.Skip);
        cmap = cmap(idx(~mem), :);
    else
        cmap = cmap(1:args.End, :);
    end
elseif isfield(args, 'Skip')
    idx = 1:size(cmap,1);
    mem = ismember(idx, args.Skip);
    cmap = cmap(idx(~mem), :);
end

if isfield(args, 'Step')
    cmap = cmap(1:args.Step:end, :);
end

if isfield(args, 'Total')
    values = round(linspace(1, size(cmap,1), args.Total));
    cmap = cmap(values, :);
end

if isfield(args, 'Repeat')
    idx = repelem(1:size(cmap,1), args.Repeat);
    cmap = cmap(idx, :);
end

cmap = cmap ./ 255;

end