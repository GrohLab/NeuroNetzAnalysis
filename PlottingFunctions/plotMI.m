function h = plotMI(MI, miNames, mName, varargin)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%% Paring inputs
p = inputParser;
% Required inputs
checkMI = @(x) all(isnumeric(x));
checkMInames = @(x) all([isstring(x), iscell(x) | isstring(x),...
    numel(x) == size(MI,2)]);
checkMName = @(x) ischar(x) | isstring(x);

% Parameters
defNBin = 32;
checkNBin = @(x) all([isnumeric(x), x > 0, numel(x) == 1]);
% Optional parameter
% Axis
defAx = gobjects(1,1);
checkAx = @(x) isa(x,'matlab.graphics.axis.Axes');
% Colormap
defCMap = lines(size(MI,2));
checkCMap = @(x) all([isnumeric(x), all(size(x) == [size(MI,2),3])]);

p.addRequired('MI',checkMI);
p.addRequired('miNames',checkMInames);
p.addRequired('mName',checkMName);
p.addParameter('numBin', defNBin, checkNBin);
p.addOptional('inTheseAxes',defAx,checkAx);
p.addOptional('colorMap',defCMap, checkCMap);


p.parse(MI, miNames, mName, varargin{:});

MI = p.Results.MI;
miNames = p.Results.miNames;
mName = p.Results.mName;
ax = p.Results.inTheseAxes;
cmap = p.Results.colorMap;
Nbin = p.Results.numBin;
%% Getting things ready
if iscell(miNames)
    miNames = string(miNames);
end
if isa(ax,'matlab.graphics.GraphicsPlaceholder')
    fig = figure('Color', 'w');
    ax = axes('Parent',fig, 'Color','none','NextPlot','add');
end
histOpts = {'BinLimits',[-1,1], 'NumBin', Nbin,...
    'Normalization', 'probability', 'EdgeColor', 'none',...
    'DisplayName'}; Nmi = size(MI,2);
fnOpts = {'UniformOutput', false};
%% Plotting MI
h = arrayfun(@(x) histogram(ax, MI(:,x), histOpts{:}, miNames(x),...
    'FaceColor', cmap(x,:)), (1:Nmi)', fnOpts{:}); h = cat(1,h{:});

end

