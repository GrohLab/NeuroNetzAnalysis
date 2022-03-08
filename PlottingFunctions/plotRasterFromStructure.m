function [axs] = plotRasterFromStructure(relSpkStruct, configStruct, varargin)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%% Auxiliary variables 
fnOpts = {'UniformOutput', false};
stOpts = {'MarkerEdgeColor','none','MarkerFaceColor','k','Marker','o'};
%% Input parsing
% Required structures checking functions
p = inputParser;
checkRSTS = @(x) isstruct(x) &...
    all(contains(fieldnames(x),{'name','SpikeTimes'}));
checkCS = @(x) isstruct(x) & all(contains(fieldnames(x),...
    {'Experiment','Viewing_window_s','Response_window_s','BinSize_s',...
    'Trigger','ConsideredConditions'}));
if checkRSTS(relSpkStruct)
    % Condition number
    Ncond = length(relSpkStruct);
    % Number of clusters and triggers (Ncl, Na)
    [Ncl, Na] = arrayfun(@(x) size(x.SpikeTimes), relSpkStruct);

    % Selecting clusters ! Assuming uniform number of units through the
    % structures
    defClSubs = (1:Ncl(1))';
    checkSubs = @(x) isPositiveIntegerValuedNumeric(x) & all(x <= Ncl(1));
    p.addParameter('ClSubs',defClSubs, checkSubs);

    % Selecting conditions
    defCond = (1:Ncond)';
    checkCond = @(x) isPositiveIntegerValuedNumeric(x) & all(x <= Ncond);
    p.addParameter('CondSubs', defCond, checkCond);
end
% 
p.addRequired('relSpkStruct', checkRSTS);
p.addRequired('configStruct', checkCS);
p.parse(relSpkStruct, configStruct, varargin{:});

relSpkStruct = p.Results.relSpkStruct;
clSubs = p.Results.ClSubs;
condSubs = p.Results.CondSubs;
configStruct = p.Results.configStruct;
%%
Nccond = length(condSubs); Nccl = length(clSubs);
mNa = min(Na); NaSubs = repmat((1:mNa)', 1, Ncond);
if any(Na > mNa)
    NaSubs = arrayfun(@(x) sort(randsample(x, mNa)), Na, fnOpts{:});
end

rastFig = figure("Name", "Raster plot", "Color", "w");
axOpts = {'NextPlot', 'Add', 'Parent', rastFig};
csp = 1; axs = gobjects(prod([Nccond, Nccl]),1);
condi = 1;
for ccond = reshape(condSubs, 1, [])
    for ccl = reshape(clSubs, 1, [])
        axs(csp) = subplot(Ncond, Nccl, csp, axOpts{:});
        XYmat = cellfun(@(x, y) [repmat(x, length(y), 1), y(:)], num2cell(1:mNa),...
            relSpkStruct(ccond).SpikeTimes(ccl, NaSubs(:,condi)), fnOpts{:});
        XYmat = cat(1, XYmat{:}); 
        scatter(axs(csp), XYmat(:,2), XYmat(:,1), stOpts{:}); 
        title(axs(csp), sprintf("clSub:%d", ccl))
        csp = csp + 1;
    end
    condi = condi + 1;
end
set(axs, "Box", "off", "Color", "none", "YAxisLocation", "origin",...
    "XLim", configStruct.Viewing_window_s);
yticks(axs, mNa); yticklabels(axs, mNa)
txOpts = {"Rotation", 90, "HorizontalAlignment","center",...
    "VerticalAlignment","bottom"}; linkaxes(axs); ylim(axs, [0,mNa+1])
arrayfun(@(x,y) text(axs(x), configStruct.Viewing_window_s(1), ...
    mNa/2, y.name, txOpts{:}), (0:Nccond-1)*Nccl + 1, relSpkStruct(condSubs))
end