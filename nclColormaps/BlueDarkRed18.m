function [cmap] = BlueDarkRed18(args)
% Edited: Nov. 11th, 2021
%
% Colormap Reference:
% https://www.ncl.ucar.edu/Document/Graphics/ColorTables/BlueDarkRed18.shtml
%
% ncolors = 18
% cmap = BlueDarkRed18()
% cmap = BlueDarkRed18(argumentName, argumentValue)
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

cmap = [[36, 0, 216];...
    [24, 28, 247];...
    [40, 87, 255];...
    [61, 135, 255];...
    [86, 176, 255];...
    [117, 211, 255];...
    [153, 234, 255];...
    [188, 249, 255];...
    [234, 255, 255];...
    [255, 255, 234];...
    [255, 241, 188];...
    [255, 214, 153];...
    [255, 172, 117];...
    [255, 120, 86];...
    [255, 61, 61];...
    [247, 39, 53];...
    [216, 21, 47];...
    [165, 0, 33]];

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