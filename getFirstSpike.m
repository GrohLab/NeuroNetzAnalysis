function [fSpkStr] = getFirstSpike(relSpkTmsStr, confStr, varargin)
%gGETFIRSTSPIKE Obvious algorithm 
%   Detailed explanation goes here
%% Auxiliary variables
fnOpts = {'UniformOutput', false};
shftOne = @(x) [x(1:end-1);x(2:end)];
foSts = @(x) [median(x), mean(x), std(x)];
%% Input parser
p = inputParser;

checkRelSpkTmsStruct = @(x) isstruct(x) & all(contains(fieldnames(x), ...
    {'name', 'SpikeTimes'}));
checkConfStr = @(x) isstruct(x) & all(contains(fieldnames(x), ...
    {'Experiment', 'Viewing_window_s', 'Response_window_s', 'BinSize_s', ...
    'Trigger', 'ConsideredConditions'}));
defRes = 5e-4;
checkRes = @(x) isnumeric(x) & numel(x) == 1 & x > 0;

if checkRelSpkTmsStruct(relSpkTmsStr)
    % Condition number
    Ncond = length(relSpkTmsStr);
    % Number of clusters and triggers (Ncl, Na)
    [Ncl, Na] = arrayfun(@(x) size(x.SpikeTimes), relSpkTmsStr);

    % Selecting clusters ! Assuming uniform number of units through the
    % structures
    defClSubs = (1:Ncl(1))';
    checkSubs = @(x) isPositiveIntegerValuedNumeric(x) & all(x <= Ncl(1));
    p.addParameter('ClSubs', defClSubs, checkSubs);

    % Selecting conditions
    defCond = (1:Ncond)';
    checkCond = @(x) isPositiveIntegerValuedNumeric(x) & all(x <= Ncond);
    p.addParameter('CondSubs', defCond, checkCond);
end
p.addRequired('relSpkTmsStr', checkRelSpkTmsStruct);
p.addRequired('confStr', checkConfStr);
p.addParameter('Res', defRes, checkRes);
p.parse(relSpkTmsStr, confStr, varargin{:});
% Results from parsing
relSpkTmsStr = p.Results.relSpkTmsStr;
confStr = p.Results.confStr;
clSubs = p.Results.ClSubs;
condSubs = p.Results.CondSubs;
res = p.Results.Res;

%% First spike for selected clusters in selected conditions
hsOpts = {'BinLimits', confStr.Response_window_s, 'BinWidth', res};
% Assuming the spike train to be a row vector, this logical operation
% returns the first spike after the trigger onset. It is not suited to find
% the first spike of a burst in a spike train
gFS = @(x) [x(~isempty(x)), xor(x(1:end-1), x(2:end))];
% Find the first spike after the trigger onset per selected condition, per
% selected cluster.
fsCell = arrayfun(@(x) cellfun(@(y) y(gFS(y>0)), x.SpikeTimes(clSubs,:), fnOpts{:}), ...
    relSpkTmsStr(condSubs), fnOpts{:});
% Trial-collapsed version for statistical purposes
fsCol = cellfun(@(x) arrayfun(@(y) cat(2, x{y,:}), (1:size(x,1))', fnOpts{:}), ...
    fsCell, fnOpts{:});
% First order statistics
% IDEA: This can be a stand alone function for other spike collection.
fsFOS = cellfun(@(x) cellfun(foSts, x, fnOpts{:}), fsCol, fnOpts{:});
fsFOS = cellfun(@(x) cat(1, x{:}), fsFOS, fnOpts{:});
% Histogram for the first spike
[fsH, hbe] = cellfun(@(x) cellfun(@(y) histcounts(y, hsOpts{:}), x, fnOpts{:}), ...
    fsCol, fnOpts{:}); hbc = mean(shftOne(hbe{1}{1}));
% Mahalanobis distance in matrix form
fsMh = cellfun(@(x) cellfun(@(y) squareform(pdist(y(:), "mahalanobis")), x, ...
    fnOpts{:}), fsCol, fnOpts{:});
end