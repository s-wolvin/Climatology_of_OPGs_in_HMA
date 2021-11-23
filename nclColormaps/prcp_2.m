function [cmap] = prcp_2(args)
% Edited: Nov. 11th, 2021
%
% Colormap Reference:
% https://www.ncl.ucar.edu/Document/Graphics/ColorTables/prcp_2.shtml
%
% ncolors = 12
% cmap = prcp_2()
% cmap = prcp_2(argumentName, argumentValue)
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

cmap = [[245, 245, 245];...
    [175, 237, 237];...
    [152, 251, 152];...
    [67, 205, 128];...
    [59, 179, 113];...
    [250, 250, 210];...
    [255, 255, 0];...
    [255, 164, 0];...
    [255, 0, 0];...
    [205, 55, 0];...
    [199, 20, 133];...
    [237, 130, 237]];

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