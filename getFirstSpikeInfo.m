function [fSpkStr] = getFirstSpikeInfo(relSpkTmsStr, confStr, varargin)
%gGETFIRSTSPIKE Obvious algorithm 
%   Detailed explanation goes here
%% Auxiliary variables
fnOpts = {'UniformOutput', false};
% shftOne = @(x) [x(1:end-1);x(2:end)];

%% Input parser
p = inputParser;

checkRelSpkTmsStruct = @(x) isstruct(x) & all(contains(fieldnames(x), ...
    {'name', 'SpikeTimes'}));

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
p.addRequired('confStr', @checkConfigStruct);
p.addParameter('Res', defRes, checkRes);
p.parse(relSpkTmsStr, confStr, varargin{:});
% Results from parsing
relSpkTmsStr = p.Results.relSpkTmsStr;
confStr = p.Results.confStr;
clSubs = p.Results.ClSubs;
condSubs = p.Results.CondSubs;
res = p.Results.Res;

%% First spike for selected clusters in selected conditions

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
[stStruct, spkDom] = getBasicSts(fsCol, confStr, Na, res);
% Getting the first spike PDF
PDF = arrayfun(@(x) cat(1, x.PDF{:}), stStruct, fnOpts{:});
fSpkStr = struct('FirstSpikeTimes', fsCell, 'FSTcollapsed', fsCol, ...
    'MedMeanStd', {stStruct.MedMeanStd}, 'PDF', PDF, 'PDF_domain', spkDom, ...
    'Quartiles', {stStruct.Quartiles}, 'MahalDist', {stStruct.MahalDist});
end
