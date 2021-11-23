function[cmap] = wind_17lev(args)
% Edited: Nov. 11th, 2021
%
% Colormap Reference:
% https://www.ncl.ucar.edu/Document/Graphics/ColorTables/wind_17lev.shtml
%
% ncolors = 18
% cmap = wind_17lev()
% cmap = wind_17lev(argumentName, argumentValue)
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

cmap = [[255, 255, 255];...
    [239, 244, 209];...
    [232, 244, 158];...
    [170, 206, 99];...
    [226, 237, 22];...
    [255, 237, 0];...
    [255, 237, 130];...
    [244, 209, 127];...
    [237, 165, 73];... 
    [229, 140, 61];...
    [219, 124, 61];...
    [239, 7, 61];...
    [232, 86, 163];...
    [155, 112, 168];...
    [99, 112, 247];...
    [127, 150, 255];...
    [142, 178, 255];... 
    [181, 201, 255]];

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