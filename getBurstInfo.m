function [outputArg1,outputArg2] = getBurstInfo(spkSubs, fs, varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Parsing input
p = inputParser;

checkFs = @(x) all([isnumeric(x), numel(x) == 1, x > 0]);

defNbin = 128;
checkNbin = @(x) all([isnumeric(x), numel(x) == 1, (round(x) - x) == 0,...
    x > 0]);

defMxDom = 1e3;
checkMxDom = @(x) all([isnumeric(x), 1/fs < x, x > 0, numel(x) == 1]);

addRequired(p, 'spkSubs', @iscell);
addRequired(p, 'fs', checkFs);
addOptional(p, 'Nbin', defNbin, checkNbin);
addOptional(p, 'MxDom', defMxDom, checkMxDom);

p.parse(spkSubs, fs, varargin{:});

spkSubs = p.Results.spkSubs;
fs = p.Results.fs;
Nbin = p.Results.Nbin;
mxDom = p.Results.MxDom;

%% Calculating ISIs and their log-distributions
fnOpts = {'UniformOutput', false};
hOpts = {'Normalization', 'probability'};
% Preparing the logarithmic domain
[isiCent, isiEdg, ~, logTau] = prepareLogBinEdges([1/fs, mxDom], Nbin);
% Getting the spike timing
spkTms = cellfun(@(x) x./fs, spkSubs, fnOpts{:});
% Getting the inter-spike interval per cluster
isiTms = cellfun(@(x) diff(x), spkTms, fnOpts{:});
% Getting the distribution per cluster in log domain
hisi = cellfun(@(x) histcounts(log10(x), isiEdg, hOpts{:}), isiTms, ...
    fnOpts{:});
hisi = cat(1,hisi{:});
% Getting the peaks and valleys of the distribution bigger than 1 ms
[hcp, xsl] = getWaveformCriticalPoints(hisi', 1/logTau);
hcp = cellfun(@(x) x - log10(fs), hcp(:,1), fnOpts{:});
msThFlags = cellfun(@(x) x > -3, hcp, fnOpts{:});
hcp = cellfun(@(x,y) x(y), hcp, msThFlags, fnOpts{:});
xsl = cellfun(@(x,y) x(y), xsl(:,1), msThFlags, fnOpts{:});
hcp = cellfun(@(x,y) x(y<0), hcp, xsl, fnOpts{:});
ycp = arrayfun(@(x) interp1(isiCent, hisi(x,:), hcp{x}), (1:size(spkSubs,1))',...
    fnOpts{:});
% Getting initial points for kmeans per cluster
[~, ordSub] = cellfun(@(x) sort(x, 'descend'), ycp, fnOpts{:});


end

